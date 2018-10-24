import os
import sys
import math
import shlex
import subprocess

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from scipy import interpolate
from scipy.integrate import quad

# ------------------------------------------------------------
# Input

PTfile = sys.argv[1]
directory = os.path.dirname(os.path.abspath(sys.argv[1]))+"/"

# Wavelength range [micron]
wavelengthMin = 0.5
wavelengthMax = 20.

# Mean molecular weight [g/mol]
mmw = 2.3

# Depolarization factor
depolarization = 0.0

# ------------------------------------------------------------

dataDir = os.path.dirname(os.path.abspath(__file__))[:-6]+'/dat/molecules/'

if not os.path.exists(directory+'opacity/'):
    os.makedirs(directory+'opacity/')

pressure, temperature = np.loadtxt(PTfile, unpack=True)

def getOpacity(filenumber):
    
    # Input:
    # filenumber corresponds with number in PTgrid.dat
    
    # Output:
    # Wavelength and opacity array
    
    wl, op = np.loadtxt(dataDir+'opacity_aver_'+str(int(filenumber)).zfill(4)+'.dat', unpack=True)
    return wl, op

def getPT(log_pressure_layer, temp_layer):
    
    # Input:
    # logarithm of the pressure in [bar]
    # temperature in [K]
    
    # Output:
    # indices array with the pressure and temperature boundaries in PTgrid.dat
    
    # Note that indices start at zero and filenames at 1
    
    opacityPTlist = np.genfromtxt(dataDir+'PTgrid.dat', skip_header=1)
    pressure_layer = 10.0**log_pressure_layer
    
    index_opac = opacityPTlist[::,0]
    p_opac = opacityPTlist[::,1]
    t_opac = opacityPTlist[::,2]
    
    for i in range(len(index_opac)):
        if t_opac[i] == temp_layer:
            upperT = t_opac[i]
            lowerT = t_opac[i]
            break
        elif t_opac[i] > temp_layer:
            upperT = t_opac[i]
            lowerT = t_opac[i-1]
            break
        elif t_opac[len(index_opac)-1] < temp_layer:
            upperT = t_opac[len(index_opac)-1]
            lowerT = t_opac[len(index_opac)-2]
            break
    
    if pressure_layer > np.max(p_opac):
        for i in range(len(index_opac)):
            if t_opac[i] == upperT:
                upperPupperTindex = i
                lowerPupperTindex = i-1
            elif t_opac[i] == lowerT:
                upperPlowerTindex = i
                lowerPlowerTindex = i-1
        indices = [upperPupperTindex, lowerPupperTindex, upperPlowerTindex, lowerPlowerTindex]
        return indices
    
    for i in range(len(index_opac)):
        if t_opac[i] == upperT:
            if p_opac[i] == pressure_layer:
                upperPupperTindex = i
                lowerPupperTindex = i
                break
            elif p_opac[i] > pressure_layer:
                upperPupperTindex = i
                lowerPupperTindex = i-1
                break
    
    for i in range(len(index_opac)):
        if t_opac[i] == lowerT:
            if p_opac[i] == pressure_layer:
                upperPlowerTindex = i
                lowerPlowerTindex = i
                break
            elif p_opac[i] > pressure_layer:
                upperPlowerTindex = i
                lowerPlowerTindex = i-1
                break
    
    try:
        indices = [upperPupperTindex, lowerPupperTindex, upperPlowerTindex, lowerPlowerTindex]
        x= p_opac[indices[0]], p_opac[indices[1]], p_opac[indices[2]], p_opac[indices[3]]
        y= t_opac[indices[0]], t_opac[indices[1]], t_opac[indices[2]], t_opac[indices[3]]
        return indices
    except UnboundLocalError:
        return [0,0,0,0]

def interpolateOpacityFiles(pressure_layer, temp_layer, indices, opacity_array):
    
    # Input:
    # pressure in [bar]
    # temperature in [K]
    # indices array of P/T boundaries in PTgrid1060.dat
    # opacity_array is a numpy array of opacities for each PT point in our 4 indexed points length is 4 x len(freqs)
    
    # Output:
    # opacity array
    
    opacityPTlist = np.genfromtxt(dataDir+'PTgrid.dat', skip_header=1)
    
    p1 = opacityPTlist[indices[1]][1]
    p2 = opacityPTlist[indices[0]][1]
    t1 = opacityPTlist[indices[2]][2]
    t2 = opacityPTlist[indices[0]][2]
    
    p1 = np.log10(p1)
    p2 = np.log10(p2)
    t1 = np.log10(t1)
    t2 = np.log10(t2)
    
    opacity_array = np.log10(opacity_array)
    temp_layer = np.log10(temp_layer)
    pressure_layer = np.log10(pressure_layer)
    
    opacity_array[opacity_array < -500] = -500
    
    if ((p1==p2) and (t1==t2)):
        return 10.**opacity_array[0]
    elif (p1==p2) :
        opacities_interp =  opacity_array[2]  + (opacity_array[0] -opacity_array[2]) * (temp_layer-t1)/(t2-t1)
        return 10.**opacities_interp
    elif (t1==t2) :
        opacities_interp =  opacity_array[1]  + (opacity_array[0] -opacity_array[1]) * (pressure_layer-p1)/(p2-p1)
        return 10.**opacities_interp
    
    R1 = ((p2-pressure_layer) / (p2-p1)) * opacity_array[3] + ((pressure_layer-p1)/(p2-p1)) * opacity_array[2]
    R2 = ((p2-pressure_layer) / (p2-p1)) * opacity_array[1] + ((pressure_layer-p1)/(p2-p1)) * opacity_array[0]
    opacities_interp = ((t2-temp_layer) / (t2-t1)) * R1 + ((temp_layer-t1)/(t2-t1)) * R2
    return 10.**opacities_interp

# Opacities

pressureLog = np.log10(pressure)

for i in range(len(pressureLog)):
    
    indices = getPT(pressureLog[i], temperature[i])

    # Just to count the number of wavelengths
    wavelength, opacity = getOpacity(indices[0]+1)
    
    opacityTotal1 = np.zeros(len(wavelength))
    opacityTotal2 = np.zeros(len(wavelength))
    opacityTotal3 = np.zeros(len(wavelength))
    opacityTotal4 = np.zeros(len(wavelength))
    opacityNew    = np.zeros(len(wavelength))
    
    wavelength, opacity1 = getOpacity(indices[0]+1)
    wavelength, opacity2 = getOpacity(indices[1]+1)
    wavelength, opacity3 = getOpacity(indices[2]+1)
    wavelength, opacity4 = getOpacity(indices[3]+1)
    
    opacity = [opacity1, opacity2, opacity3, opacity4]
    opacityNew = interpolateOpacityFiles(pressure[i], temperature[i], indices, opacity)
    
    # Individual molecule opacity and total absorption for each layer
    absFile = directory+"opacity/absorption_"+str('%02d' % (len(pressureLog)-i))+".dat"
    f = open(absFile,'w')
    f.write("# Wavelength [micron] - Opacity x VMR [cm2/molecule]\n\n")
    for i in range(len(wavelength)):
        f.write(str(wavelength[i])+'\t'+str(opacityNew[i])+'\n')
    f.close()

# Make FITS files

opacityDir = directory+'opacity/'

avogadro = 6.02214129e23
loschmidt = 2.6867805e19 # [cm-3]

def refractiveIndexH2(wavelength):
    
	a = 13.58e-5
	b = 7.52e-3
	ri = 1. + a + a*b / (wavelength*wavelength)
	
	return ri

def rayleighP11(theta):
    
	alpha = math.cos(theta)
	delta = (1.-depolarization) / (1.+depolarization/2.)
	ray = ( ((alpha*alpha + 1.) * delta) +  (1.-delta) ) * math.sin(theta)
	
	return ray

def rayleighScatter(alpha):
    
	delta      = (1.-depolarization) / (1.+depolarization/2.)
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

layer = 0

for file in os.listdir(opacityDir):
    if file.startswith('absorption_') and file.endswith('.dat'):
    	
		wavelength, absorption = np.loadtxt(opacityDir+file, unpack=True)
		mass = mmw / avogadro # Molecule mass [g]
		absorption /= mass # [cm2 molecule-1] -> [cm2 g-1]

		scatter = np.zeros((180,16,len(wavelength)))

		opacityWavelength = []
		opacityScatter = []
		opacityAbsorption = []
		opacityExtinction = []

		wl = 0.

		def rayleighH2(wavelength):
			sigma = (8.14e-13/(wavelength**4)) + (1.28e-6/(wavelength**6)) + (1.61/(wavelength**8))
			return sigma

		for i in range(len(wavelength)):
			if wavelength[i] >= wavelengthMin:

				# Rayleigh cross section [cm2]
				ri = refractiveIndexH2(wavelength[i])
				rindex = (ri*ri-1.)*(ri*ri-1.) / ((ri*ri+2.)*(ri*ri+2.))
				dep = ( 6.+3.*depolarization ) / ( 6.-7.*depolarization )
				crossSection = 24.*math.pi*math.pi*math.pi*rindex*dep / ( ((wavelength[i]*1.e-4)**4) * (loschmidt**2) )

				opacityWavelength.append(wavelength[i])
				opacityScatter.append(crossSection/mass)
				opacityAbsorption.append(absorption[i])
				opacityExtinction.append(crossSection/mass+absorption[i])

				wl = wavelength[i]

				if wavelength[i] > wavelengthMax:
					break

		opacity = np.zeros((4,len(opacityWavelength)))
		scatter = np.zeros((180,16,len(opacityWavelength)))

		for i in range(len(opacityWavelength)):
			opacity[0,i] = opacityWavelength[i] # [micron]
			opacity[1,i] = opacityExtinction[i] # [cm2 g-1]
			opacity[2,i] = opacityAbsorption[i] # [cm2 g-1]
			opacity[3,i] = opacityScatter[i] # [cm2 g-1]

		# Scattering matrices

		scatter = np.zeros((180,16,len(opacityWavelength)))

		for i in range(len(opacityWavelength)):
			for j in range(180):

				scatterMatrixLow = rayleighScatter(math.cos(float(j)*math.pi/180.))
				scatterMatrixUp = rayleighScatter(math.cos(float(j+1)*math.pi/180.))

				for m in range(16):

					scatter[j,m,i] = ( scatterMatrixLow[m] + scatterMatrixUp[m] ) / 2.
					scatter[j,m,i] /= rayleighNorm

		fitsfile = opacityDir+'gas_opacity_'+file[11:13]+'.fits'
		
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

for file in os.listdir(opacityDir):
	if file.startswith('absorption_') and file.endswith('.dat'):
		os.remove(opacityDir+file)
