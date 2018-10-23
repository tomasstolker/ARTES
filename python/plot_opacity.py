import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import sys, os
from astropy.io import fits

# ------------------------------------------------------------
# Input

atmosphere = sys.argv[1]

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))
opacityDir = directory[:-6]+'input/'+atmosphere+'/opacity/'
plotDir = directory[:-6]+'input/'+atmosphere+'/plot/'

count = 0
    
for file in os.listdir(opacityDir):
    if file.endswith('.fits') and not "gas_opacity_" in file[:-5]:
        count = count + 1

colors = iter(plt.cm.rainbow(np.linspace(0,1,count)))

j = 0

for file in os.listdir(opacityDir):
    if file.endswith('.fits'):

        fitsfile = opacityDir+file
        
        hdulist = fits.open(fitsfile)
        hdu0 = hdulist[0]
        hdu1 = hdulist[1]
        hdulist.close()

        naxis = hdu0.header['NAXIS1']

        wavelength = np.zeros(naxis)
        absorption = np.zeros(naxis)
        scattering = np.zeros(naxis)

        data = fits.getdata(fitsfile,0)

        for i in range(naxis):
            wavelength[i] = data[0,i]
            absorption[i] = data[2,i]
            scattering[i] = data[3,i]

        if not "gas_opacity_" in file[:-5]:

            c = next(colors)

            if np.size(wavelength) == 1:
                plt.plot(wavelength, scattering, color=c, marker='o', ms=5, mew=0, label=file[:-5])
                plt.plot(wavelength, absorption, color=c, marker='^', ms=5, mew=0)
            else:
                plt.plot(wavelength, scattering, color=c, ls='--', label=file[:-5])
                plt.plot(wavelength, absorption, color=c, ls='-')
                
                j += 1

        else:

            if np.size(wavelength) == 1:
                plt.plot(wavelength, scattering, color='gray', marker='o', ms=3, mew=0)
                plt.plot(wavelength, absorption, color='gray', marker='^', ms=3, mew=0)
            else:
                plt.plot(wavelength, scattering, color='gray', ls='--', lw=0.6)
                plt.plot(wavelength, absorption, color='gray', ls='-', lw=0.6)

plt.xlabel('Wavelength [micron]')
plt.ylabel('Opacity [cm$^2$/g]')

plt.xscale('log')
plt.yscale('log')

if count != j:
    plt.legend(loc='upper left')

plt.savefig(plotDir+'opacity.pdf', bbox_inches='tight')