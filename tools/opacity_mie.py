import numpy as np
import os, shutil, math, sys
from astropy.io import fits
from scipy.integrate import simps
from sys import platform as _platform

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

fitsfile = 'forsterite.fits'
riFile = 'forsterite.dat'

nr = 1000 # Number of radius points
nf = 20 # Number of volume fractions for distribution of hollow spheres (DHS)
density = 1. # Particle density [g cm-3]
amin = 0.01 # Minimum particle size [micron]
amax = 2.0 # Maximum particle size [micron]
apow = 0.0 # Size distribution power law index
fmax = 0.0 # Irregularity parameter, maximum volume void fraction for DHS, fmax=0. -> Mie
r_eff = 0.1 # Effective radius [micron], overrules amin, amax, apow
v_eff = 0.05 # Effective variance [dimensionless], overrules amin, amax, apow

# Wavelengths from refractive index file
riWavelength = False

# Wavelengths from FITS file
fitsWavelength = False
wavelengthFile = ''

# Wavelength range [micron]
manualWavelength = True
wavelengthMin = 0.7
wavelengthMax = 0.7
step = 0.2

# ------------------------------------------------------------

scriptDir =  os.path.dirname(os.path.abspath(__file__))
riDir = scriptDir[:-6]+'dat/refractive/'
fitsfile = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+fitsfile
wavelengthFile = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+wavelengthFile
tempDir = scriptDir[:-6]+'temp/'

# Percentage always set to 100%, i.e. no mixture of multiple grain species
percentage = 100.

if _platform == "linux" or _platform == "linux2":
    ComputePart = scriptDir[:-6]+'bin/ComputePartLinux'
elif _platform == "darwin":
    ComputePart = scriptDir[:-6]+'bin/ComputePartMac'
elif _platform == "win32":
    sys.exit("Error: ComputePartWindows not available yet")
else:
    sys.exit("Error: Operating system not recognized, ", _platform)

if not os.path.exists(tempDir):
    os.makedirs(tempDir)

os.chdir(tempDir)
	
wavelength, real, imaginary = np.loadtxt(riDir+riFile, unpack=True)

f = open("wavelength.dat",'w')

if riWavelength:

	if np.size(wavelength) == 1:
		f.write(str(wavelength))

	elif np.size(wavelength) > 1:

		for i in range(len(wavelength)):
			f.write(str(wavelength[i])+"\n")
	
elif fitsWavelength:
	
	hdulist = fits.open(wavelengthFile)
	hdu = hdulist[0].data
	wavelength = hdu[0]
	hdulist.close()
	for i in range(len(wavelength)):
		f.write(str(wavelength[i])+"\n")
	
elif manualWavelength:

	for i in range(int((wavelengthMax-wavelengthMin)/step)+1):
		f.write(str(wavelengthMin+i*step)+"\n")
	
f.close()

f = open("mie.in",'w')
f.write(str(nr)+"\n")
f.write(str(nf)+"\n")
f.write("'"+str(riDir+riFile)+"' \n")
f.write(str(percentage)+"\t"+str(density)+"\t"+str(amin)+"\t"+ \
	str(amax)+"\t"+str(apow)+"\t"+str(fmax))
f.close()

os.chmod(ComputePart, 0700)
if r_eff > 0.:
    print ComputePart+" mie.in wavelength.dat "+str(r_eff)+" "+str(v_eff)
    os.system(ComputePart+" mie.in wavelength.dat "+str(r_eff)+" "+str(v_eff))
else:
	os.system(ComputePart+" mie.in wavelength.dat")
os.rename(tempDir+"particle.fits", fitsfile)
shutil.rmtree(tempDir)

hdulist = fits.open(fitsfile)
wavelengths = int(hdulist[1].header['NAXIS1'])
elements = int(hdulist[1].header['NAXIS2'])
opacity = hdulist[0].data
scatter = hdulist[1].data
hdulist.close()

# Change to 16-elements scatter matrix
	
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

# Scatter matrix normalization

angle = np.zeros(180)
for i in range(180):
	angle[i] = (float(i)+0.5) * math.pi/180.

for j in range(wavelengths):

	p11Int = scatter[:,0,j]
	norm = simps(p11Int*np.sin(angle), angle)
	norm *= 2.*math.pi
	scatter[:,:,j] /= norm

# Write updated FITS file

hdulist = fits.HDUList()
hdulist.append(fits.ImageHDU(np.array(opacity), name='opacity'))
hdulist.append(fits.ImageHDU(np.array(scatter), name='scatter'))
hdu = hdulist[0]
hdu.header['COMMENT'] = '1. Wavelength [micron]'
hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
hdulist.writeto(fitsfile, overwrite=True)
hdulist.close()