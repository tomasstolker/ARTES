import matplotlib.pyplot as plt
import numpy as np
import sys, os
from astropy.io import fits

# ------------------------------------------------------------
# Input

output = sys.argv[1]
azimuth = int(sys.argv[2])
scale = float(sys.argv[3])

# ------------------------------------------------------------

directory =  os.path.dirname(os.path.abspath(__file__))
plotDir = directory[:-6]+'output/'+output+'/plot/'

hdulist = fits.open(directory[:-6]+'output/'+output+'/input/atmosphere.fits')
r_face = hdulist[0].data # [m]
t_face = hdulist[1].data # [deg]
p_face = hdulist[2].data # [deg]
hdulist.close()

r_face *= 1e-3 # [km]
r_face -= r_face[0]
t_face = 90. - t_face

r = np.zeros((np.size(r_face)-1))
t = np.zeros((np.size(t_face)-1))

for i in range(np.size(r_face)-1):
    r[i] = (r_face[i]+r_face[i+1])/2.

for i in range(np.size(t_face)-1):
    t[i] = (t_face[i]+t_face[i+1])/2.

tt, rr = np.meshgrid(t,r)

hdulist = fits.open(directory[:-6]+'output/'+output+'/output/global.fits')
data = hdulist[0].data[azimuth]
flow_r = data[:,:,0]
flow_t = data[:,:,1]
flow_p = data[:,:,2]
hdulist.close()

flow_r = np.transpose(flow_r)
flow_t = np.transpose(flow_t)
flow_p = np.transpose(flow_p)

fig = plt.quiver(tt, rr, flow_t, flow_r, flow_p, pivot='middle', cmap='inferno', clim=[-1,1], scale=scale)

cb = plt.colorbar(fig)
cb.set_label('Azimuthal direction')

plt.xlim(90,-90)
plt.ylim(0.,max(r_face))

plt.xticks(np.arange(90,-100,-30))

plt.xlabel('Latitude [deg]')
plt.ylabel('Height [km]')

plt.savefig(plotDir+'global.pdf', bbox_inches='tight')