import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import sys, os
from astropy.io import fits
from matplotlib.colors import LogNorm

# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")
azimuth = int(sys.argv[2])

# ------------------------------------------------------------

wavelength, cell_depth = np.loadtxt(output+'output/cell_depth.dat', unpack=True)

cell_depth = int(cell_depth)

# Read grid cell faces

hdulist = fits.open(output+'input/atmosphere.fits')
r_face = hdulist[0].data # [m]
t_face = hdulist[1].data # [deg]
p_face = hdulist[2].data # [deg]
hdulist.close()

r_face *= 1e-3 # [km]
r_face -= r_face[0]
t_face = 90. - t_face # [deg]

# Pressure-latitude meshgrid

tt,rr = np.meshgrid(t_face, r_face[cell_depth+1:])

# Read Latitudinal flow

hdulist = fits.open(output+'output/latitudinal.fits')
data = hdulist[0].data[azimuth]
flow_up = data[:,:,0]
flow_down = data[:,:,1]
flow_south = data[:,:,2]
flow_north = data[:,:,3]
hdulist.close()

# Transpose

flow_up = np.transpose(flow_up)
flow_down = np.transpose(flow_down)
flow_south = np.transpose(flow_south)
flow_north = np.transpose(flow_north)

# Cell depth and higher

flow_up = flow_up[cell_depth:-1,:]
flow_down = flow_down[cell_depth+1:,:]
flow_south = flow_south[cell_depth+1:,:]
flow_north = flow_north[cell_depth+1:,:]

# Minimum and maximum

NS_min = min( np.min(flow_south[np.nonzero(flow_south)]), np.min(flow_north[np.nonzero(flow_north)]) )
NS_max = max( np.max(flow_south[np.nonzero(flow_south)]), np.max(flow_north[np.nonzero(flow_north)]) )

# Plot up

plt.figure(figsize=(8,5))

fig = plt.pcolormesh(tt, rr, flow_up, cmap='inferno', norm=LogNorm(), edgecolor='black', linewidth=0.5)

cb = plt.colorbar(fig)
cb.ax.set_ylabel(r'L$_{\rm up}$/L$_{\rm atm}$')

plt.xlim(90,-90)
plt.ylim(np.amin(rr),np.amax(rr))

plt.xticks(np.arange(90,-100,-30))

plt.xlabel('Latitude [deg]')
plt.ylabel('Altitude [km]')

plt.savefig(output+'/plot/flow_up.pdf', bbox_inches='tight')
plt.clf()

# Plot down

plt.figure(figsize=(8,5))

fig = plt.pcolormesh(tt, rr, flow_down, cmap='inferno', norm=LogNorm(), edgecolor='black', linewidth=0.5)

cb = plt.colorbar(fig)
cb.ax.set_ylabel(r'L$_{\rm down}$/L$_{\rm atm}$')

plt.xlim(90,-90)
plt.ylim(np.amin(rr),np.amax(rr))

plt.xticks(np.arange(90,-100,-30))

plt.xlabel('Latitude [deg]')
plt.ylabel('Altitude [km]')

plt.savefig(output+'/plot/flow_down.pdf', bbox_inches='tight')
plt.clf()

# Plot south

plt.figure(figsize=(8,5))

fig = plt.pcolormesh(tt, rr, flow_south, cmap='inferno', vmin=NS_min, vmax=NS_max, norm=LogNorm(), edgecolor='black', linewidth=0.5)

cb = plt.colorbar(fig)
cb.ax.set_ylabel(r'L$_{\rm south}$/L$_{\rm atm}$')

plt.xlim(90,-90)
plt.ylim(np.amin(rr),np.amax(rr))

plt.xticks(np.arange(90,-100,-30))

plt.xlabel('Latitude [deg]')
plt.ylabel('Altitude [km]')

plt.savefig(output+'/plot/flow_south.pdf', bbox_inches='tight')
plt.clf()

# Plot north

plt.figure(figsize=(8,5))

fig = plt.pcolormesh(tt, rr, flow_north, cmap='inferno', vmin=NS_min, vmax=NS_max, norm=LogNorm(), edgecolor='black', linewidth=0.5)

cb = plt.colorbar(fig)
cb.ax.set_ylabel(r'L$_{\rm north}$/L$_{\rm atm}$')

plt.xlim(90,-90)
plt.ylim(np.amin(rr),np.amax(rr))

plt.xticks(np.arange(90,-100,-30))

plt.xlabel('Latitude [deg]')
plt.ylabel('Altitude [km]')

plt.savefig(output+'/plot/flow_north.pdf', bbox_inches='tight')
plt.clf()

# Plot ratio

nr = np.size(flow_north,0)
ntheta = np.size(flow_north,1)

diff = np.zeros((nr,ntheta))

for i in range(nr):
    for j in range(ntheta):
        
        if j == 0:
            diff[i,j] = flow_south[i,j]-flow_north[i,j+1]
        elif j > 0 and j < ntheta-1:
            diff[i,j] = (flow_north[i,j]+flow_south[i,j])-(flow_north[i,j+1]+flow_south[i,j-1])
        elif j == ntheta-1:
            diff[i,j] = flow_north[i,j]-flow_south[i,j-1]

        # diff[i,j] /= flow_up[i,j]-flow_down[i,j]

cmap = {name:plt.get_cmap(name) for name in ('copper', 'bone_r')}
N = 50
vmin = np.amin(diff)
vmax = np.amax(diff)
levels = np.concatenate([np.linspace(vmin, 0, N, endpoint=False), np.linspace(0, vmax, N+1, endpoint=True)])
colors = np.concatenate([cmap[name](np.linspace(0, 1, N)) for name in ('bone_r', 'copper')])
cmap, norm = mcolors.from_levels_and_colors(levels, colors)

plt.figure(figsize=(8,5))

fig = plt.pcolormesh(tt, rr, diff, cmap=cmap, norm=norm, edgecolor='black', linewidth=0.5)

cb = plt.colorbar(fig)
cb.ax.set_ylabel(r'$\Delta$L$_{\rm hor}$/$\Delta$L$_{\rm ver}$')

plt.xlim(90,-90)
plt.ylim(np.amin(rr),np.amax(rr))

plt.xticks(np.arange(90,-100,-30))

plt.xlabel('Latitude [deg]')
plt.ylabel('Altitude [km]')

plt.savefig(output+'/plot/flow_ratio.pdf', bbox_inches='tight')
