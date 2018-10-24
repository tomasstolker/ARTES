import os
import sys
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from astropy.io import fits

# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")
scaling = float(sys.argv[2])
labelPol = float(sys.argv[3]) # [%]

# ------------------------------------------------------------

fitsDir = output+'output/'
plotDir = output+'plot/'

hdulist = fits.open(fitsDir+'stokes.fits')
data = hdulist[0].data
stokesI = data[0,:,:]
stokesQ = data[1,:,:]
stokesU = data[2,:,:]
stokesV = data[3,:,:]
hdulist.close()

hdulist = fits.open(fitsDir+'error.fits')
data = hdulist[0].data
errorI = data[0,:,:]
errorQ = data[1,:,:]
errorU = data[2,:,:]
errorV = data[3,:,:]
errorP = data[4,:,:]
hdulist.close()

npix = np.size(stokesI,0)

polarization = np.zeros((npix,npix))
for i in range(npix):
    for j in range(npix):
        if stokesI[i,j] > 0.:
            polarization[i,j] = np.sqrt(stokesQ[i,j]**2+stokesU[i,j]**2)/stokesI[i,j]

direction = 0.5*np.arctan2(stokesU, stokesQ)

# Stokes I

fig = plt.imshow(stokesI, origin='lower', cmap='gist_heat', interpolation='nearest', vmin=0., vmax=np.amax(stokesI), extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plotDir+'image_I.pdf', bbox_inches='tight')
plt.clf()

# Stokes Q

Qmax = np.amax(np.abs(stokesQ))
Qmin = -Qmax

cmap = {name:plt.get_cmap(name) for name in ('copper', 'bone_r')}
N = 50
levels = np.concatenate([np.linspace(Qmin, 0, N, endpoint=False), np.linspace(0, Qmax, N+1, endpoint=True)])
colors = np.concatenate([cmap[name](np.linspace(0, 1, N)) for name in ('bone_r', 'copper')])
cmap, norm = mcolors.from_levels_and_colors(levels, colors)

fig = plt.imshow(stokesQ, origin='lower', cmap=cmap, norm=norm, interpolation='nearest', vmin=Qmin, vmax=Qmax, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plotDir+'image_Q.pdf', bbox_inches='tight')
plt.clf()

# Stokes U

Umax = np.amax(np.abs(stokesU))
Umin = -Umax

cmap = {name:plt.get_cmap(name) for name in ('copper', 'bone_r')}
N = 50
levels = np.concatenate([np.linspace(Umin, 0, N, endpoint=False), np.linspace(0, Umax, N+1, endpoint=True)])
colors = np.concatenate([cmap[name](np.linspace(0, 1, N)) for name in ('bone_r', 'copper')])
cmap, norm = mcolors.from_levels_and_colors(levels, colors)

fig = plt.imshow(stokesU, origin='lower', cmap=cmap, norm=norm, interpolation='nearest', vmin=Umin, vmax=Umax, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plotDir+'image_U.pdf', bbox_inches='tight')
plt.clf()

# Stokes V

if abs(np.amax(stokesV)) > 0.:

    Vmax = np.amax(stokesV)
    Vmin = np.amin(stokesV)

    fig = plt.imshow(stokesV, origin='lower', cmap='gist_heat', interpolation='nearest', vmin=Vmin, vmax=Vmax, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

    cb = plt.colorbar(fig)
    cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

    plt.xlabel('x offset', fontsize=12)
    plt.ylabel('y offset', fontsize=12)

    plt.savefig(plotDir+'image_V.pdf', bbox_inches='tight')
    plt.clf()

# Error Stokes I

errorI_frac = np.zeros((npix,npix))
for i in range(npix):
    for j in range(npix):
        if polarization[i,j] > 0.:
            errorI_frac[i,j] = errorI[i,j]/stokesI[i,j]

errorI_max = np.nanmax(errorI_frac)
if errorI_max > 1.:
    errorI_max = 1.

fig = plt.imshow(errorI_frac, origin='lower', cmap='magma', interpolation='nearest', vmin=0.0, vmax=errorI_max, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Fractional error Stokes I')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plotDir+'error_I.pdf', bbox_inches='tight')
plt.clf()

# Error degree of polarization

errorP_frac = np.zeros((npix,npix))
for i in range(npix):
    for j in range(npix):
        if polarization[i,j] > 0.:
            errorP_frac[i,j] = errorP[i,j]/polarization[i,j]

errorP_max = np.nanmax(errorP_frac)
if errorP_max > 1.:
    errorP_max = 1.

fig = plt.imshow(errorP_frac, origin='lower', cmap='magma', interpolation='nearest', vmin=0.0, vmax=errorP_max, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Fractional error degree of polarization')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plotDir+'error_P.pdf', bbox_inches='tight')
plt.clf()

# Polarization map

fig = plt.imshow(stokesI, origin='lower', cmap='gist_heat', interpolation='nearest', vmin=0., vmax=np.amax(stokesI), extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

k = int(npix/50.)+1
m = int(k/2)

for i in range(npix):
    for j in range(npix):
        if polarization[j,i] > 0. and i%k == m and j%k == m:
            plt.quiver(float(i)-npix/2.+0.5, float(j)-npix/2.+0.5,
                       math.cos(direction[j,i]+math.pi), math.sin(direction[j,i]+math.pi),
                       angles='xy', scale=scaling/polarization[j,i], color='#66CCFF',
                       edgecolor='#66CCFF', pivot='middle', headwidth=0, headlength=0, headaxislength=0, width=0.002)

label = str("{:.0f}".format(labelPol))+'%'
plt.quiver(npix/15.-npix/2., npix/20.-npix/2., 1., 0., scale=scaling/(labelPol/100.), width=0.007, pivot='middle', headwidth=0, headlength=0, headaxislength=0, color='#66CCFF', edgecolor='#66CCFF')
plt.annotate(label, xy=(npix/15.-npix/2., npix/20.-npix/2.), xycoords='data', color='white', xytext=(0,7), textcoords='offset points', ha='center')

deg = np.sqrt(np.sum(stokesQ)**2+np.sum(stokesU)**2)/np.sum(stokesI)
polDeg = str("{:.2f}".format(100.*deg))+'%'
plt.annotate(polDeg, xy=(0,0), xycoords='axes fraction', xytext=(0.97,0.97), textcoords='axes fraction', color='white', weight='bold', ha='right', va='top')

plt.savefig(plotDir+'polarization.pdf', bbox_inches='tight')