import numpy as np
import math, os, sys
from astropy.io import fits
from scipy.integrate import quad

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

fitsOutput = 'rayleigh.fits'

MMWscat = 2.02 # Mean molecular weight of scattering gas [g/mol]
depolarization = 0.0
singleScatteringAlbedo = 1.0

# Wavelengths from FITS file
fitsWavelength = False
wavelengthFile = ''

# Wavelength range [micron]
manualWavelength = True
wavelengthMin = 0.7
wavelengthMax = 0.7
step = 1.0

# ------------------------------------------------------------

scriptDir =  os.path.dirname(os.path.abspath(__file__))
fitsOutput = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+fitsOutput
wavelengthFile = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+wavelengthFile

if fitsWavelength:
    
    hdulist = fits.open(wavelengthFile)
    hdu = hdulist[0].data
    wavelength = hdu[0]
    hdulist.close()

elif manualWavelength:

    wavelength = []
    for i in range(int((wavelengthMax-wavelengthMin)/step)+1):
        wavelength.append(wavelengthMin+float(i)*step)

avogadro = 6.02214129e23
loschmidt = 2.6867805e19 # [cm-3]
gasMassScat = MMWscat / avogadro # Molecule mass [g]

opacityWavelength = []
opacityScatter = []
opacityAbsorption = []
opacityExtinction = []

for i in range(len(wavelength)):
    
    # Rayleigh cross section [cm2]

    a = 13.58e-5
    b = 7.52e-3
    ri = 1. + a + a*b / (wavelength[i]*wavelength[i]) # Wavelength in [micron]

    rindex = (ri*ri-1.)*(ri*ri-1.) / ((ri*ri+2.)*(ri*ri+2.))
    dep = ( 6.+3.*depolarization ) / ( 6.-7.*depolarization )
    crossSection = 24.*math.pi*math.pi*math.pi*rindex*dep / ( ((wavelength[i]*1.e-4)**4) * (loschmidt**2) )

    opacityScat = crossSection/gasMassScat

    opacityWavelength.append(wavelength[i])
    opacityScatter.append(opacityScat)
    opacityExtinction.append(opacityScat/singleScatteringAlbedo)
    opacityAbsorption.append(opacityScat/singleScatteringAlbedo-opacityScat)

opacity = np.zeros((4,len(opacityWavelength)))

for i in range(len(opacityWavelength)):
    opacity[0,i] = opacityWavelength[i]
    opacity[1,i] = opacityExtinction[i]
    opacity[2,i] = opacityAbsorption[i]
    opacity[3,i] = opacityScatter[i]

def rayleighP11(theta):

    alpha = math.cos(theta)
    delta = (1.-depolarization) / (1.+depolarization/2.)
    ray = ( ((alpha*alpha + 1.) * delta) +  (1.-delta) ) * math.sin(theta)

    return ray

rayleighNorm, error = quad(rayleighP11, 0., math.pi)
rayleighNorm *= 2.*math.pi

def rayleighScatter(alpha):

    scatterMatrix = np.zeros((16))
    
    delta = (1.-depolarization) / (1.+depolarization/2.)
    deltaPrime = (1.-2.*depolarization) / (1.-depolarization)

    scatterMatrix[0] = alpha*alpha + 1.
    scatterMatrix[1] = alpha*alpha - 1.
    scatterMatrix[4] = scatterMatrix[1]
    scatterMatrix[5] = scatterMatrix[0]
    scatterMatrix[10] = 2.*alpha
    scatterMatrix[15] = deltaPrime*scatterMatrix[10]

    scatterMatrix = delta * scatterMatrix
    scatterMatrix[0] = scatterMatrix[0] + (1.-delta)

    return scatterMatrix

scatter = np.zeros((180,16,len(opacityWavelength)))

for i in range(len(opacityWavelength)):
    for j in range(180):

        scatterMatrixLow = rayleighScatter(math.cos(float(j)*math.pi/180.))
        scatterMatrixUp = rayleighScatter(math.cos(float(j+1)*math.pi/180.))

        for m in range(16):

            scatter[j,m,i] = ( scatterMatrixLow[m] + scatterMatrixUp[m] ) / 2.
            scatter[j,m,i] /= rayleighNorm

hdulist = fits.HDUList()
hdulist.append(fits.ImageHDU(np.array(opacity), name='opacity'))
hdulist.append(fits.ImageHDU(np.array(scatter), name='scatter'))
hdu = hdulist[0]
hdu.header['COMMENT'] = '1. Wavelength [micron]'
hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
hdulist.writeto(fitsOutput, overwrite=True)
hdulist.close()