import os
import sys
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from astropy.io import fits


# ------------------------------------------------------------
# Input

output = sys.argv[1]
scaling = float(sys.argv[2])
labelPol = float(sys.argv[3])  # [%]

# ------------------------------------------------------------

fits_dir = os.path.join(output, 'output/')
plot_dir = os.path.join(output, 'plot/')

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

data = fits.getdata(fits_dir+'stokes.fits')

stokes_i = data[0, ]
stokes_q = data[1, ]
stokes_u = data[2, ]
stokes_v = data[3, ]

data = fits.getdata(fits_dir+'error.fits')
error_i = data[0, ]
error_q = data[1, ]
error_u = data[2, ]
error_v = data[3, ]
error_p = data[4, ]

npix = stokes_i.shape[0]

polarization = np.zeros((npix, npix))
for i in range(npix):
    for j in range(npix):
        if stokes_i[i,j] > 0.:
            polarization[i,j] = np.sqrt(stokes_q[i, j]**2+stokes_u[i, j]**2)/stokes_i[i, j]

direction = 0.5*np.arctan2(stokes_u, stokes_q)

# Stokes I

fig = plt.imshow(stokes_i, origin='lower', cmap='gist_heat', interpolation='nearest', vmin=0., vmax=np.amax(stokes_i), extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plot_dir+'image_i.pdf', bbox_inches='tight')
plt.clf()

# Stokes Q

Qmax = np.amax(np.abs(stokes_q))
Qmin = -Qmax

cmap = {name:plt.get_cmap(name) for name in ('copper', 'bone_r')}
N = 50
levels = np.concatenate([np.linspace(Qmin, 0, N, endpoint=False), np.linspace(0, Qmax, N+1, endpoint=True)])
colors = np.concatenate([cmap[name](np.linspace(0, 1, N)) for name in ('bone_r', 'copper')])
cmap, norm = mcolors.from_levels_and_colors(levels, colors)

fig = plt.imshow(stokes_q, origin='lower', cmap=cmap, norm=norm, interpolation='nearest', vmin=Qmin, vmax=Qmax, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plot_dir+'image_q.pdf', bbox_inches='tight')
plt.clf()

# Stokes U

Umax = np.amax(np.abs(stokes_u))
Umin = -Umax

cmap = {name:plt.get_cmap(name) for name in ('copper', 'bone_r')}
N = 50
levels = np.concatenate([np.linspace(Umin, 0, N, endpoint=False), np.linspace(0, Umax, N+1, endpoint=True)])
colors = np.concatenate([cmap[name](np.linspace(0, 1, N)) for name in ('bone_r', 'copper')])
cmap, norm = mcolors.from_levels_and_colors(levels, colors)

fig = plt.imshow(stokes_u, origin='lower', cmap=cmap, norm=norm, interpolation='nearest', vmin=Umin, vmax=Umax, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plot_dir+'image_u.pdf', bbox_inches='tight')
plt.clf()

# Stokes V

if abs(np.amax(stokes_v)) > 0.:

    Vmax = np.amax(stokes_v)
    Vmin = np.amin(stokes_v)

    fig = plt.imshow(stokes_v, origin='lower', cmap='gist_heat', interpolation='nearest', vmin=Vmin, vmax=Vmax, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

    cb = plt.colorbar(fig)
    cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

    plt.xlabel('x offset', fontsize=12)
    plt.ylabel('y offset', fontsize=12)

    plt.savefig(plot_dir+'image_v.pdf', bbox_inches='tight')
    plt.clf()

# Error Stokes I

error_i_frac = np.zeros((npix, npix))
for i in range(npix):
    for j in range(npix):
        if polarization[i,j] > 0.:
            error_i_frac[i,j] = error_i[i,j]/stokes_i[i,j]

error_i_max = np.nanmax(error_i_frac)
if error_i_max > 1.:
    error_i_max = 1.

fig = plt.imshow(error_i_frac, origin='lower', cmap='magma', interpolation='nearest', vmin=0.0, vmax=error_i_max, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Fractional error Stokes I')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plot_dir+'error_i.pdf', bbox_inches='tight')
plt.clf()

# Error degree of polarization

error_p_frac = np.zeros((npix, npix))
for i in range(npix):
    for j in range(npix):
        if polarization[i,j] > 0.:
            error_p_frac[i,j] = error_p[i, j]/polarization[i, j]

error_p_max = np.nanmax(error_p_frac)
if error_p_max > 1.:
    error_p_max = 1.

fig = plt.imshow(error_p_frac, origin='lower', cmap='magma', interpolation='nearest', vmin=0.0, vmax=error_p_max, extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Fractional error degree of polarization')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

plt.savefig(plot_dir+'error_p.pdf', bbox_inches='tight')
plt.clf()

# Polarization map

fig = plt.imshow(stokes_i, origin='lower', cmap='gist_heat', interpolation='nearest', vmin=0., vmax=np.amax(stokes_i), extent=[-npix/2.,npix/2.,-npix/2.,npix/2.])

cb = plt.colorbar(fig)
cb.ax.set_ylabel('Surface brightness [W m$^{-2}$ micron$^{-1}$ mas$^{-2}$]')

plt.xlabel('x offset', fontsize=12)
plt.ylabel('y offset', fontsize=12)

k = int(npix/50.)+1
m = int(k/2)

for i in range(npix):
    for j in range(npix):
        if polarization[j, i] > 0. and error_p_frac[j, i] < 0.1 and i%k == m and j%k == m:
            plt.quiver(float(i)-npix/2.+0.5, float(j)-npix/2.+0.5,
                       math.cos(direction[j,i]+math.pi), math.sin(direction[j,i]+math.pi),
                       angles='xy', scale=scaling/polarization[j,i], color='#66CCFF',
                       edgecolor='#66CCFF', pivot='middle', headwidth=0, headlength=0, headaxislength=0, width=0.002)

label = str('{:.0f}'.format(labelPol))+'%'
plt.quiver(npix/15.-npix/2., npix/20.-npix/2., 1., 0., scale=scaling/(labelPol/100.), width=0.007, pivot='middle', headwidth=0, headlength=0, headaxislength=0, color='#66CCFF', edgecolor='#66CCFF')
plt.annotate(label, xy=(npix/15.-npix/2., npix/20.-npix/2.), xycoords='data', color='white', xytext=(0,7), textcoords='offset points', ha='center')

deg = np.sqrt(np.sum(stokes_q)**2+np.sum(stokes_u)**2)/np.sum(stokes_i)
polDeg = str('{:.2f}'.format(100.*deg))+'%'
plt.annotate(polDeg, xy=(0,0), xycoords='axes fraction', xytext=(0.97,0.97), textcoords='axes fraction', color='white', weight='bold', ha='right', va='top')

plt.savefig(plot_dir+'polarization.pdf', bbox_inches='tight')
