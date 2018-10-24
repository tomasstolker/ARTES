import numpy as np
import math, os, sys
from astropy.io import fits
from scipy.integrate import quad

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

fitsOutput = 'gas.fits'
absorptionCoefficients = 'methane.dat'

VMR = 1.8e-3 # Volume mixing ratio of absorbing molecule
MMWabs = 16.04 # Mean molecular weight of absorbing gas [g/mol]
MMWscat = 2.02 # Mean molecular weight of scattering gas [g/mol]
depolarization = 0.02

# Wavelength range [micron]
wavelengthMin = 0.4
wavelengthMax = 1.0

# Wavelengths from absorption coefficients
absorptionWavelength = False

# Manual wavelengths
manualWavelength = True
manualStep = 0.01

# ------------------------------------------------------------

scriptDir =  os.path.dirname(os.path.abspath(__file__))
fitsOutput = scriptDir[:-6]+'input/'+atmosphere+'/opacity/'+fitsOutput
absorptionCoefficients = scriptDir[:-6]+'dat/absorption/'+absorptionCoefficients

avogadro = 6.02214129e23 # [mol-1]
loschmidt = 2.6867805e19 # [cm-3]

gasMassAbs = MMWabs / avogadro # Molecule mass [g]
gasMassScat = MMWscat / avogadro # Molecule mass [g]

# Wavelength [micron] - Absorption [cm2 mol-1]
w, a = np.loadtxt(absorptionCoefficients, unpack=True)

a = a / gasMassAbs # [cm2 molecule-1] -> [cm2 g-1]

opacityWavelength = []
opacityScatter = []
opacityAbsorption = []
opacityExtinction = []

wavelength = []
absorption = []

if absorptionWavelength:

    for i in range(len(w)):

        if w[i] >= wavelengthMin:
            wavelength.append(w[i])
            absorption.append(a[i])
            
        if w[i] > wavelengthMax:
            break

elif manualWavelength:

    wl = wavelengthMin

    for i in range(len(w)):

        if w[i] >= wl and w[i] < wl+manualStep:
            
            wavelength.append(w[i])
            absorption.append(a[i])
            wl += manualStep
            
        if w[i] > wavelengthMax:
            break

for i in range(len(wavelength)):

    # H2 refractive index

    a = 13.58e-5
    b = 7.52e-3 # [micron^2]
    ri = 1. + a + a*b / (wavelength[i]*wavelength[i])

    # Rayleigh cross section [cm2]	

    dep = ( 6.+3.*depolarization ) / ( 6.-7.*depolarization )
    rindex = ((ri*ri-1.)/loschmidt)**2.
    crossSection = (8.*math.pi**3./3.)*rindex*dep
    crossSection /= (wavelength[i]*1.e-4)**4.

    opacityWavelength.append(wavelength[i])
    opacityScatter.append(crossSection/gasMassScat)
    opacityAbsorption.append(absorption[i]*VMR)
    opacityExtinction.append(crossSection/gasMassScat+absorption[i]*VMR)

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

def rayleighScatter(alpha):
	
	delta = (1.-depolarization) / (1.+depolarization/2.)
	deltaPrime = (1.-2.*depolarization) / (1.-depolarization)

	scatterMatrix     = np.zeros((16))
	scatterMatrix[0]  = alpha*alpha + 1.
	scatterMatrix[1]  = alpha*alpha - 1.
	scatterMatrix[4]  = scatterMatrix[1]
	scatterMatrix[5]  = scatterMatrix[0]
	scatterMatrix[10] = 2.*alpha
	scatterMatrix[15] = deltaPrime*scatterMatrix[10]
	scatterMatrix     = delta * scatterMatrix
	scatterMatrix[0]  = scatterMatrix[0] + (1.-delta)

	return scatterMatrix

rayleighNorm, error = quad(rayleighP11, 0., math.pi)
rayleighNorm *= 2.*math.pi

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