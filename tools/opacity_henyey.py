import os
import sys
import math

import numpy as np

from astropy.io import fits
from scipy.integrate import quad

# ------------------------------------------------------------
# Input

opacity = sys.argv[1]

fitsOutput = 'henyey.fits'

absorption = 0.05 # Absorption opacity [cm2 g-1]
scattering = 1.0 # Scattering opacity [cm2 g-1]

# Henyey-Greenstein parameters
g1 = 0.8
w1 = 1.0
g2 = 0.
w2 = 0.
g3 = 0.
w3 = 0.
pLinear = 0.5
pCircular = 0.0
skew = 0.0

# Wavelengths from FITS file
fitsWavelength = True
wavelengthFile = 'ammonia_a1.0_b0.10.fits'

# Wavelength range [micron]
manualWavelength = False
wavelengthMin = 0.7
wavelengthMax = 0.7
step = 1.0

# ------------------------------------------------------------

fitsOutput = os.path.join(opacity, fitsOutput)
wavelengthFile = os.path.join(opacity, wavelengthFile)

if fitsWavelength:
    
    hdulist = fits.open(wavelengthFile)
    hdu = hdulist[0].data
    wavelength = hdu[0]
    hdulist.close()

elif manualWavelength:

    wavelength = []
    for i in range(int((wavelengthMax-wavelengthMin)/step)+1):
        wavelength.append(wavelengthMin+float(i)*step)

opacityScatter = []
opacityExtinction = []

opacity = np.zeros((4,len(wavelength)))
for i in range(len(wavelength)):
    opacity[0,i] = wavelength[i]
    opacity[1,i] = absorption + scattering
    opacity[2,i] = absorption
    opacity[3,i] = scattering

def hgP11(theta):
    alpha = math.cos(theta)
    henyey = w1 * (1.-g1*g1) / ((1.+g1*g1-2.*g1*alpha)**(3./2.))
    henyey += w2 * (1.-g2*g2) / ((1.+g2*g2-2.*g2*alpha)**(3./2.))
    henyey += w3 * (1.-g3*g3) / ((1.+g3*g3-2.*g3*alpha)**(3./2.))
    henyey *= math.sin(theta)

    return henyey

def HGscatter(alpha):
    scatterMatrix = np.zeros((16))
    
    alphaF = alpha * (1.+ 3.13 * skew * math.exp(-7.*alpha/math.pi))
    cosAlphaF = math.cos(alphaF)

    scatterMatrix[0] = w1 * (1.-g1*g1) / ((1.+g1*g1-2.*g1*alpha)**(3./2.))
    scatterMatrix[0] += w2 * (1.-g2*g2) / ((1.+g2*g2-2.*g2*alpha)**(3./2.))
    scatterMatrix[0] += w3 * (1.-g3*g3) / ((1.+g3*g3-2.*g3*alpha)**(3./2.))
    scatterMatrix[1] = -pLinear * scatterMatrix[0] * (1.-alpha*alpha) / (1.+alpha*alpha)
    scatterMatrix[4] = scatterMatrix[1]
    scatterMatrix[5] = scatterMatrix[0]
    scatterMatrix[10] = scatterMatrix[0] * (2.*alpha) / (1.+alpha*alpha)
    scatterMatrix[11] = pCircular * scatterMatrix[5] * (1.-cosAlphaF*cosAlphaF) / (1.+cosAlphaF*cosAlphaF)
    scatterMatrix[14] = -scatterMatrix[11]
    scatterMatrix[15] = scatterMatrix[10]

    return scatterMatrix

hgNorm, error = quad(hgP11, 0., math.pi)
hgNorm *= 2.*math.pi

scatter = np.zeros((180, 16, len(wavelength)))

for i in range(len(wavelength)):
    for j in range(180):
        scatterMatrixLow = HGscatter(math.cos(float(j)*math.pi/180.))
        scatterMatrixUp = HGscatter(math.cos(float(j+1)*math.pi/180.))

        for m in range(16):
            scatter[j, m, i] = (scatterMatrixLow[m] + scatterMatrixUp[m]) / 2.
            scatter[j, m, i] /= hgNorm

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
