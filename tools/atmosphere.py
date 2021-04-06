import configparser
import math
import os
import sys

import numpy as np

from astropy.io import fits
from scipy.integrate import simps


config_file = os.path.abspath(sys.argv[1])
directory = os.path.dirname(config_file) + '/'

# Normalize phase functions

opacity_dir = directory+'opacity/'

angle = np.zeros(180)
for i in range(180):
    angle[i] = (float(i)+0.5) * math.pi/180.

for file in os.listdir(opacity_dir):
    if file.endswith('.fits'):

        fitsfile = opacity_dir+file
        hdulist = fits.open(fitsfile)
        wavelengths = int(hdulist[1].header['NAXIS1'])
        elements = int(hdulist[1].header['NAXIS2'])
        opacity = hdulist[0].data
        scatter = hdulist[1].data
        hdulist.close()

        # Change 6-elements to 16-elements scattering matrix

        if elements == 6:

            scatter_new = np.zeros((180, 16, wavelengths))

            for i in range(wavelengths):
                for j in range(180):

                    scatter_new[j, 0, i] = scatter[j, 0, i]
                    scatter_new[j, 1, i] = scatter[j, 1, i]
                    scatter_new[j, 4, i] = scatter[j, 1, i]
                    scatter_new[j, 5, i] = scatter[j, 2, i]
                    scatter_new[j, 10, i] = scatter[j, 3, i]
                    scatter_new[j, 11, i] = scatter[j, 4, i]
                    scatter_new[j, 14, i] = -scatter[j, 4, i]
                    scatter_new[j, 15, i] = scatter[j, 5, i]

            scatter = scatter_new
    
        for j in range(wavelengths):

            p11_int = scatter[:,0,j]
            norm = simps(p11_int*np.sin(angle), angle)
            norm *= 2.*math.pi
            scatter[:, :, j] /= norm

        # Write updated FITS file

        try:
            os.remove(fitsfile)
        except OSError:
            pass

        hdulist = fits.HDUList()
        hdulist.append(fits.ImageHDU(np.array(opacity), name='opacity'))
        hdulist.append(fits.ImageHDU(np.array(scatter), name='scattermatrix'))
        hdu = hdulist[0]
        hdu.header['COMMENT'] = '1. Wavelength [um]'
        hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
        hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
        hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
        hdulist.writeto(fitsfile)
        hdulist.close()

# Read input file

parser = configparser.ConfigParser()
parser.read(config_file)

# Pressure temperature profile

radial = []
theta = []
phi = []

gas_constant = 8.3144621 # [J K-1 mol-1]

# Planet radius [Rjup]
r_planet = float(parser.get('grid', 'radius'))
r_planet *= 69911e3 # [m]

# Ring
ring = parser.has_option('composition', 'ring')

# Gas opacities [on/off]
gas = parser.getboolean('composition', 'gas')

density_gas = []

if os.path.isfile(directory+'pressure_temperature.dat'):

    # Mean molecular mass [kg/mol]
    mmw = float(parser.get('composition', 'molweight')) * 1.e-3 # [g/mol] -> [kg/mol]
    # Surface gravity
    logg = float(parser.get('composition', 'log_g'))
    gravity = 1e-2*10.**logg # [m s-2]

    # Pressure temperature profile
    pressure, temperature = np.loadtxt(directory+'pressure_temperature.dat', unpack=True)

    if isinstance(pressure, np.float32) or isinstance(pressure, np.float64):
        pressure = np.array([pressure])
        temperature = np.array([temperature])

    # [bar] -> [Pa]
    pressure *= 1.e5
    # Reverse order
    pressure = pressure[::-1]
    temperature = temperature[::-1]

    scale_height = np.zeros(len(pressure))
    density_gas = np.zeros(len(pressure))
    radial = np.zeros(len(pressure))

    radial[0] = 0.
    scale_height[0] =  gas_constant * temperature[0] / ( mmw * gravity )
    density_gas[0] = pressure[0] / ( gravity * scale_height[0] )

    if len(pressure) == 1:
        raise ValueError('Minimum 2 pressures are required in pressure_temperature.dat')

    for i in range(1, len(pressure)):

        scale_height[i] =  gas_constant * temperature[i] / ( mmw * gravity ) # [m]
        density_gas[i] = pressure[i] / ( gravity * scale_height[i] ) # [kg m-3]
        radial[i] = radial[i-1] - scale_height[i] * np.log( pressure[i] / pressure[i-1] ) # [m]

        if radial[i]-radial[i-1] < 0.:
            print('Radial error')

    # Number of radial cell faces, starting to count at 1
    nr = len(radial)

    pressure = pressure[:-1]
    temperature = temperature[:-1]
    scale_height = scale_height[:-1]
    density_gas = density_gas[:-1]
    radialGas = radial[:-1]
    
else:
    
    r = parser.get('grid', 'radial')
    rr = [ chunk.strip() for chunk in r.split(',') ]
    
    # Number of radial cell faces, starting to count at 1
    nr = len(rr)+1
    
    radial.append(0.)
    if nr > 0 and len(rr[0]) > 0:
        for i in range(int(nr-1)):
            radial.append(float(rr[i])*1.e3) # [km] -> [m]

for i in range(len(radial)):
    radial[i] += r_planet

# Polar and azimuthal cell faces

t = parser.get('grid', 'theta')
tt = [ chunk.strip() for chunk in t.split(',') ]

# Number of latitudinal cell faces, starting to count at 1
if not tt[0]:
    ntheta = 2
else:
    ntheta = int(len(tt)+2)

theta.append(0.)
if len(tt) > 0 and tt[0]:
    for i in range(len(tt)):
        theta.append(float(tt[i]))
theta.append(180.)

p = parser.get('grid', 'phi')
pp = [ chunk.strip() for chunk in p.split(',') ]

# Number of azimuthal cell faces, starting to count at 1
if not pp[0]:
    nphi = 1
else:
    nphi = int(len(pp)+1)

phi.append(0.)
if len(pp) > 0 and pp[0]:
    for i in range(len(pp)):
        phi.append(float(pp[i]))

# Read gas opacity files

if gas:

    opacity = fits.getdata(directory+'opacity/gas_opacity_01.fits')
    nwav = np.size(opacity, 1)
    wavelengths = opacity[0]

    opacity_gas = np.zeros((np.size(density_gas), 4, nwav))
    scatter_gas = np.zeros((np.size(density_gas), 180, 16, nwav))
    
    for i in range(np.size(density_gas)):
        gas_fits = directory+'opacity/gas_opacity_'+str(i+1).zfill(2)+'.fits'
        hdulist = fits.open(gas_fits)
        opacity = hdulist[0].data
        scatter = hdulist[1].data
        hdulist.close()
        
        opacity_gas[i, :, :] = opacity / 10. # [m2 kg-1]
        scatter_gas[i, :, :, :] = scatter

# Read other opacity files

opacity_num = 1
while True:

    fits_num = 'fits'+str(opacity_num).zfill(2)

    try:
        parser.get('composition', fits_num)
    except:
        configparser.NoOptionError
        break

    opacity_num += 1

opacity_num -= 1

if opacity_num > 0:
    
    fits_file = parser.get('composition', 'fits01')
    opacity = fits.getdata(directory+'opacity/'+fits_file)
    nwav = np.size(opacity, 1)
    wavelengths = opacity[0]

    opacity_other = np.zeros((opacity_num, 4, nwav))
    scatter_other = np.zeros((opacity_num, 180, 16, nwav))

    for i in range(opacity_num):

        fits_num = 'fits'+str(i+1).zfill(2)
        fits_file = directory+'opacity/'+parser.get('composition', fits_num)

        hdulist = fits.open(fits_file)
        opacity = hdulist[0].data
        scatter = hdulist[1].data
        hdulist.close()

        # Also wavelength divided by 10, but not needed
        opacity_other[i, :, :] = opacity / 10. # [m2 kg-1]
        scatter_other[i, :, :, :] = scatter

# Read composition from configuration file

composition = []
r_in = []
r_out = []
theta_in = []
theta_out = []
phi_in = []
phi_out = []
density_other = []

i = 1

while True:

    opacity_num = 'opacity'+str(i).zfill(2)

    try:
        parser.get('composition', opacity_num)
    except:
        configparser.NoOptionError
        break

    a = parser.get('composition', opacity_num)
    aa = [ chunk.strip() for chunk in a.split(',') ]

    if 'nr' in aa[3]:
        aa[3] = nr-1

    if 'ntheta' in aa[5]:
        aa[5] = ntheta-1

    if 'nphi' in aa[7]:
        aa[7] = nphi
        
    composition.append(int(aa[0]))
    r_in.append(int(aa[2]))
    r_out.append(int(aa[3]))
    theta_in.append(int(aa[4]))
    theta_out.append(int(aa[5]))
    phi_in.append(int(aa[6]))
    phi_out.append(int(aa[7]))
	
    try:
        float(aa[1])
        check = True
    except ValueError:
        check = False

    if check:
        density_other.append(float(aa[1])*1e3) # [kg m-3]

    i += 1

# Opacity x density [m-1]

opacity_scattering = np.zeros((nwav, nphi, ntheta-1, nr-1))
opacity_absorption = np.zeros((nwav, nphi, ntheta-1, nr-1))
scatter = np.zeros((180, 16, nwav, nphi, ntheta-1, nr-1))
density = np.zeros((nphi, ntheta-1, nr-1))

if gas:

    for i in range(nr-1):
        for j in range(ntheta-1):
            for k in range(nphi):
                for m in range(nwav):
                    opacity_absorption[m, k, j, i] = density_gas[i] * opacity_gas[i, 2, m]                    
                    opacity_scattering[m, k, j, i] = density_gas[i] * opacity_gas[i, 3, m]
                    scatter[:, :, m, k, j, i] = scatter_gas[i, :, :, m]

    for i in range(nr-1):
        density[:, :, i] = density_gas[i]  # [kg m-3]

if len(composition) > 0:
    
    for n in range(len(composition)):
        for m in range(nwav):
            for k in range(phi_in[n], phi_out[n]):
                for j in range(theta_in[n], theta_out[n]):
                    for i in range(r_in[n], r_out[n]):
                        
                        oScat = density_other[n]*opacity_other[composition[n]-1, 3, m]
                        oAbs = density_other[n]*opacity_other[composition[n]-1, 2, m]
                        
                        if density[k, j, i] == 0.:
                            scatter[:, :, m, k, j, i] = scatter_other[composition[n]-1, :, :, m]

                        elif density[k, j, i] > 0.:
                            weight = (oScat+oAbs)/(oScat+oAbs+opacity_scattering[m, k, j, i]+opacity_absorption[m, k, j, i])
                            scatter[:, :, m, k, j, i] *= (1.-weight)
                            scatter[:, :, m, k, j, i] += weight*scatter_other[composition[n]-1, :, :, m]

                        opacity_scattering[m, k, j, i] += oScat
                        opacity_absorption[m, k, j, i] += oAbs

    for n in range(len(composition)):
        for k in range(phi_in[n], phi_out[n]):
            for j in range(theta_in[n], theta_out[n]):
                for i in range(r_in[n], r_out[n]):
                    density[k, j, i] += density_other[composition[n]-1] # [kg m-3]


# Temperature [K]

temperature_grid = np.zeros((nphi, ntheta-1, nr-1))

if os.path.isfile(directory+'pressure_temperature.dat'):

    for k in range(nphi):
        for j in range(ntheta-1):
            for i in range(nr-1):
                temperature_grid[k, j, i] = temperature[i]

# Ring system

if ring:
    
    a = parser.get('composition', 'ring')
    aa = [ chunk.strip() for chunk in a.split(',') ]

    r_max = np.amax(radial)
    
    ring_in = r_max+float(aa[3])*1e3 # [m]
    ring_out = r_max+float(aa[4])*1e3 # [m]
    
    radial = np.append(radial, ring_in)
    radial = np.append(radial, ring_out)

    ring_density = np.zeros((nphi, ntheta-1, 2))
    ring_density[:, int(aa[5]):int(aa[6]), 1] = float(aa[1])
    density = np.append(density, ring_density, axis=2)
    
    ring_temperature = np.zeros((nphi, ntheta-1, 2))
    ring_temperature[:, int(aa[5]):int(aa[6]), 1] = float(aa[2])
    temperature_grid = np.append(temperature_grid, ring_temperature, axis=2)

    dust2gas = float(aa[7])
    gas_abs = float(aa[8]) # [cm2 g-1]
    gas_scat = 0. # [cm2 g-1]
    
    ring_opacity_scat = np.zeros((nwav, nphi, ntheta-1,  2))
    ring_opacity_abs = np.zeros((nwav, nphi, ntheta-1, 2))
    ring_scatter = np.zeros((180, 16, nwav, nphi, ntheta-1, 2))
    
    for m in range(nwav):
        oScat = dust2gas*float(aa[1])*opacity_other[int(aa[0])-1, 3, m] + (1.-dust2gas)*float(aa[1])*gas_scat
        oAbs = dust2gas*float(aa[1])*opacity_other[int(aa[0])-1, 2, m] + (1.-dust2gas)*float(aa[1])*gas_abs

        ring_opacity_scat[m, :, int(aa[5]):int(aa[6]), 1] = oScat
        ring_opacity_abs[m, :, int(aa[5]):int(aa[6]), 1] = oAbs
        
        for i in range(180):
            for j in range(16):
                ring_scatter[i, j, m, :, int(aa[5]):int(aa[6]), 1] = scatter_other[int(aa[0])-1, i, j, m]
        
    opacity_scattering = np.append(opacity_scattering, ring_opacity_scat, axis=3)
    opacity_absorption = np.append(opacity_absorption, ring_opacity_abs, axis=3)
    scatter = np.append(scatter, ring_scatter, axis=5)

# Asymmetry parameter

asymmetry = np.zeros((nwav, nphi, ntheta-1, nr-1))
for m in range(nwav):
    for k in range(nphi):
        for j in range(ntheta-1):
            for i in range(nr-1):
                g = 2.*math.pi*simps(scatter[:, 0, m, k, j, i]*np.cos(angle)*np.sin(angle), angle)
                asymmetry[m, k, j, i] = g

# Write FITS output

with fits.HDUList() as hdunew:
    hdunew.append(fits.ImageHDU(radial, name='radial'))  # [m]
    hdunew.append(fits.ImageHDU(theta, name='polar'))  # [deg]
    hdunew.append(fits.ImageHDU(phi, name='azimuthal'))  # [deg]
    hdunew.append(fits.ImageHDU(wavelengths, name='wavelength'))  # [um]
    hdunew.append(fits.ImageHDU(density, name='density'))  # [kg m-3]
    hdunew.append(fits.ImageHDU(temperature_grid, name='temperature'))  # [K]
    hdunew.append(fits.ImageHDU(opacity_scattering, name='scattering'))  # [m-1]
    hdunew.append(fits.ImageHDU(opacity_absorption, name='absorption'))  # [m-1]
    hdunew.append(fits.ImageHDU(scatter, name='scattermatrix'))
    hdunew.append(fits.ImageHDU(asymmetry, name='asymmetry'))
    hdunew.writeto(directory+'atmosphere.fits', overwrite=True)
