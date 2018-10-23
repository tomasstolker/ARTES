import numpy as np
import ConfigParser, sys, os, math
from astropy.io import fits
from scipy.integrate import simps

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

# ------------------------------------------------------------

directory = os.path.dirname(os.path.abspath(__file__))
directory = directory[:-6]+'input/'+atmosphere+'/'

if not os.path.exists(directory):
    print 'The '+atmosphere+' folder is missing!'
    answer = raw_input('Create folder '+atmosphere+'? [y/n]')
    if answer == 'y':
        if not os.path.exists(directory):
            print 'Creating '+directory
            os.makedirs(directory)
        if not os.path.exists(directory+'opacity'):
            print 'Creating '+directory+'opacity'
            os.makedirs(directory+'opacity')
        if not os.path.exists(directory+'plot'):
            print 'Creating '+directory+'plot'
            os.makedirs(directory+'plot')
    sys.exit(0)

try:
    os.remove(directory+'atmosphere.fits')
except OSError:
    pass

# Normalize phase functions

opacityDir = directory+'opacity/'

angle = np.zeros(180)
for i in range(180):
    angle[i] = (float(i)+0.5) * math.pi/180.

for file in os.listdir(opacityDir):
    if file.endswith(".fits"):

        fitsfile = opacityDir+file
        hdulist = fits.open(fitsfile)
        wavelengths = int(hdulist[1].header['NAXIS1'])
        elements = int(hdulist[1].header['NAXIS2'])
        opacity = hdulist[0].data
        scatter = hdulist[1].data
        hdulist.close()

        # Change 6-elements to 16-elements scattering matrix

        if elements == 6:

            scatterNew = np.zeros((180,16,wavelengths))

            for i in range(wavelengths):
                for j in range(180):

                    scatterNew[j,0,i] = scatter[j,0,i]
                    scatterNew[j,1,i] = scatter[j,1,i]
                    scatterNew[j,4,i] = scatter[j,1,i]
                    scatterNew[j,5,i] = scatter[j,2,i]
                    scatterNew[j,10,i] = scatter[j,3,i]
                    scatterNew[j,11,i] = scatter[j,4,i]
                    scatterNew[j,14,i] = -scatter[j,4,i]
                    scatterNew[j,15,i] = scatter[j,5,i]

            scatter = scatterNew
    
        for j in range(wavelengths):

            p11Int = scatter[:,0,j]
            norm = simps(p11Int*np.sin(angle), angle)
            norm *= 2.*math.pi
            scatter[:,:,j] /= norm

        # Write updated FITS file

        try:
            os.remove(fitsfile)
        except OSError:
            pass

        hdulist = fits.HDUList()
        hdulist.append(fits.ImageHDU(np.array(opacity), name='opacity'))
        hdulist.append(fits.ImageHDU(np.array(scatter), name='scattermatrix'))
        hdu = hdulist[0]
        hdu.header['COMMENT'] = '1. Wavelength [micron]'
        hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
        hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
        hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
        hdulist.writeto(fitsfile)
        hdulist.close()

# Read input file

if not os.path.isfile(directory+'atmosphere.in'):
    sys.exit('The atmosphere.in file is missing!')

parser = ConfigParser.SafeConfigParser()
parser.read(directory+'atmosphere.in')

# Pressure temperature profile

radial = []
theta = []
phi = []

gasConstant = 8.3144621 # [J K-1 mol-1]

# Planet radius [Rjup]
r_planet = float(parser.get('grid', 'radius'))
r_planet *= 69911e3 # [m]

# Ring
ring = parser.has_option('composition', 'ring')

# Gas opacities [on/off]
gas = parser.getboolean("composition", "gas")

densityGas = []

if os.path.isfile(directory+'pressure_temperature.dat'):

    # Mean molecular mass [kg/mol]
    mmw = float(parser.get('composition', 'molweight')) * 1.e-3 # [g/mol] -> [kg/mol]
    # Surface gravity
    logg = float(parser.get('composition', 'logg'))
    gravity = 1e-2*10.**logg # [m s-2]

    # Pressure temperature profile
    pressure, temperature = np.loadtxt(directory+'pressure_temperature.dat', unpack=True)
    # [bar] -> [Pa]
    pressure *= 1.e5
    # Reverse order
    pressure = pressure[::-1]
    temperature = temperature[::-1]

    scaleHeight = np.zeros(len(pressure))
    densityGas = np.zeros(len(pressure))
    radial = np.zeros(len(pressure))

    radial[0] = 0.
    scaleHeight[0] =  gasConstant * temperature[0] / ( mmw * gravity )
    densityGas[0] = pressure[0] / ( gravity * scaleHeight[0] )

    for i in range(1,len(pressure)):

        scaleHeight[i] =  gasConstant * temperature[i] / ( mmw * gravity ) # [m]
        densityGas[i] = pressure[i] / ( gravity * scaleHeight[i] ) # [kg m-3]
        radial[i] = radial[i-1] - scaleHeight[i] * np.log( pressure[i] / pressure[i-1] ) # [m]

        if radial[i]-radial[i-1] < 0.:
            print "Radial error"

    # Number of radial cell faces, starting to count at 1
    nr = len(radial)

    pressure = pressure[:-1]
    temperature = temperature[:-1]
    scaleHeight = scaleHeight[:-1]
    densityGas = densityGas[:-1]
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
    nwav = np.size(opacity,1)
    wavelengths = opacity[0]

    opacityGas = np.zeros((np.size(densityGas),4,nwav))
    scatterGas = np.zeros((np.size(densityGas),180,16,nwav))
    
    for i in range(np.size(densityGas)):
        
        gasFITS = directory+'opacity/gas_opacity_'+str(i+1).zfill(2)+'.fits'
        hdulist = fits.open(gasFITS)
        opacity = hdulist[0].data
        scatter = hdulist[1].data
        hdulist.close()
        
        opacityGas[i,:,:] = opacity / 10. # [m2 kg-1]
        scatterGas[i,:,:,:] = scatter

# Read other opacity files

opacityNumber = 1
while True:

    fitsNumber = 'fits'+str(opacityNumber).zfill(2)
    try:
        parser.get('composition', fitsNumber)
    except:
        ConfigParser.NoOptionError
        break

    opacityNumber += 1

opacityNumber -= 1

if opacityNumber > 0:
    
    fitsFile = parser.get('composition', 'fits01')
    opacity = fits.getdata(directory+'opacity/'+fitsFile)
    nwav = np.size(opacity,1)
    wavelengths = opacity[0]

    opacityOther = np.zeros((opacityNumber,4,nwav))
    scatterOther = np.zeros((opacityNumber,180,16,nwav))

    for i in range(opacityNumber):

        fitsNumber = 'fits'+str(i+1).zfill(2)
        fitsFile = directory+'opacity/'+parser.get('composition', fitsNumber)

        hdulist = fits.open(fitsFile)
        opacity = hdulist[0].data
        scatter = hdulist[1].data
        hdulist.close()

        # Also wavelength divided by 10, but not needed
        opacityOther[i,:,:] = opacity / 10. # [m2 kg-1]
        scatterOther[i,:,:,:] = scatter

# Read composition from atmosphere.in

composition = []
rIn = []
rOut = []
thetaIn = []
thetaOut = []
phiIn = []
phiOut = []
densityOther = []

i = 1
while True:
	
    opacityNumber = 'opacity'+str(i).zfill(2)
    try:
        parser.get('composition', opacityNumber)
    except:
        ConfigParser.NoOptionError
        break

    a = parser.get('composition', opacityNumber)
    aa = [ chunk.strip() for chunk in a.split(',') ]

    if 'nr' in aa[3]:
        aa[3] = nr-1
    if 'ntheta' in aa[5]:
        aa[5] = ntheta-1
    if 'nphi' in aa[7]:
        aa[7] = nphi
        
    composition.append(int(aa[0]))
    rIn.append(int(aa[2]))
    rOut.append(int(aa[3]))
    thetaIn.append(int(aa[4]))
    thetaOut.append(int(aa[5]))
    phiIn.append(int(aa[6]))
    phiOut.append(int(aa[7]))
	
    try:
        float(aa[1])
        check = True
    except ValueError:
        check = False

    if check:
        densityOther.append(float(aa[1])*1e3) # [kg m-3]

    i += 1

# Opacity x density [m-1]

opacityScattering = np.zeros((nwav, nphi, ntheta-1, nr-1))
opacityAbsorption = np.zeros((nwav, nphi, ntheta-1, nr-1))
scatter = np.zeros((180, 16, nwav, nphi, ntheta-1, nr-1))
density = np.zeros((nphi, ntheta-1, nr-1))

if gas:

    for i in range(nr-1):
        for j in range(ntheta-1):
            for k in range(nphi):
                for m in range(nwav):
                    
                    opacityAbsorption[m,k,j,i] = densityGas[i] * opacityGas[i,2,m]                    
                    opacityScattering[m,k,j,i] = densityGas[i] * opacityGas[i,3,m]
                    scatter[:,:,m,k,j,i] = scatterGas[i,:,:,m]

    for i in range(nr-1):
        density[:,:,i] = densityGas[i] # [kg m-3]

if len(composition) > 0:
    
    for n in range(len(composition)):
        for m in range(nwav):
            for k in range(phiIn[n],phiOut[n]):
                for j in range(thetaIn[n],thetaOut[n]):
                    for i in range(rIn[n],rOut[n]):
                        
                        oScat = densityOther[n]*opacityOther[composition[n]-1,3,m]
                        oAbs = densityOther[n]*opacityOther[composition[n]-1,2,m]
                        
                        if density[k,j,i] == 0.:

                            scatter[:,:,m,k,j,i] = scatterOther[composition[n]-1,:,:,m]

                        elif density[k,j,i] > 0.:

                            weight = (oScat+oAbs)/(oScat+oAbs+opacityScattering[m,k,j,i]+opacityAbsorption[m,k,j,i])
                            
                            scatter[:,:,m,k,j,i] *= (1.-weight)
                            scatter[:,:,m,k,j,i] += weight*scatterOther[composition[n]-1,:,:,m]

                        opacityScattering[m,k,j,i] += oScat
                        opacityAbsorption[m,k,j,i] += oAbs
                        
    for n in range(len(composition)):
        for k in range(phiIn[n],phiOut[n]):
            for j in range(thetaIn[n],thetaOut[n]):
                for i in range(rIn[n],rOut[n]):

                    density[k,j,i] += densityOther[composition[n]-1] # [kg m-3]


# Temperature [K]

temperatureGrid = np.zeros((nphi, ntheta-1, nr-1))

if os.path.isfile(directory+'pressure_temperature.dat'):

    for k in range(nphi):
        for j in range(ntheta-1):
            for i in range(nr-1):

                temperatureGrid[k,j,i] = temperature[i]

# Ring system

if ring:
    
    a = parser.get('composition', 'ring')
    aa = [ chunk.strip() for chunk in a.split(',') ]

    r_max = np.amax(radial)
    
    ring_in = r_max+float(aa[3])*1e3 # [m]
    ring_out = r_max+float(aa[4])*1e3 # [m]
    
    radial = np.append(radial,ring_in)
    radial = np.append(radial,ring_out)

    ring_density = np.zeros((nphi,ntheta-1,2))
    ring_density[:,int(aa[5]):int(aa[6]),1] = float(aa[1])
    density = np.append(density,ring_density,axis=2)
    
    ring_temperature = np.zeros((nphi,ntheta-1,2))
    ring_temperature[:,int(aa[5]):int(aa[6]),1] = float(aa[2])
    temperatureGrid = np.append(temperatureGrid,ring_temperature,axis=2)

    dust2gas = float(aa[7])
    gas_abs = float(aa[8]) # [cm2 g-1]
    gas_scat = 0. # [cm2 g-1]
    
    ring_opacity_scat = np.zeros((nwav,nphi,ntheta-1,2))
    ring_opacity_abs = np.zeros((nwav,nphi,ntheta-1,2))
    ring_scatter = np.zeros((180,16,nwav,nphi,ntheta-1,2))
    
    for m in range(nwav):
    
        oScat = dust2gas*float(aa[1])*opacityOther[int(aa[0])-1,3,m] + (1.-dust2gas)*float(aa[1])*gas_scat
        oAbs = dust2gas*float(aa[1])*opacityOther[int(aa[0])-1,2,m] + (1.-dust2gas)*float(aa[1])*gas_abs
        
        ring_opacity_scat[m,:,int(aa[5]):int(aa[6]),1] = oScat
        ring_opacity_abs[m,:,int(aa[5]):int(aa[6]),1] = oAbs
        
        for i in range(180):
            for j in range(16):
                ring_scatter[i,j,m,:,int(aa[5]):int(aa[6]),1] = scatterOther[int(aa[0])-1,i,j,m]
        
    opacityScattering = np.append(opacityScattering,ring_opacity_scat,axis=3)
    opacityAbsorption = np.append(opacityAbsorption,ring_opacity_abs,axis=3)
    scatter = np.append(scatter,ring_scatter,axis=5)

# Asymmetry parameter

asymmetry = np.zeros((nwav, nphi, ntheta-1, nr-1))
for m in range(nwav):
    for k in range(nphi):
        for j in range(ntheta-1):
            for i in range(nr-1):
                g = 2.*math.pi*simps(scatter[:,0,m,k,j,i]*np.cos(angle)*np.sin(angle), angle)
                asymmetry[m,k,j,i] = g

# Write FITS output

hdunew = fits.HDUList()
hdunew.append(fits.ImageHDU(radial, name='radial')) # [m]
hdunew.append(fits.ImageHDU(theta, name='polar')) # [deg]
hdunew.append(fits.ImageHDU(phi, name='azimuthal')) # [deg]
hdunew.append(fits.ImageHDU(wavelengths, name='wavelength')) # [micron]
hdunew.append(fits.ImageHDU(density, name='density')) # [kg m-3]
hdunew.append(fits.ImageHDU(temperatureGrid, name='temperature')) # [K]
hdunew.append(fits.ImageHDU(opacityScattering, name='scattering')) # [m-1]
hdunew.append(fits.ImageHDU(opacityAbsorption, name='absorption')) # [m-1]
hdunew.append(fits.ImageHDU(scatter, name='scattermatrix'))
hdunew.append(fits.ImageHDU(asymmetry, name='asymmetry'))
hdunew.writeto(directory+'atmosphere.fits', overwrite=True)
hdunew.close()