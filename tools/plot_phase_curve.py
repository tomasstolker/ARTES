import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits


# ------------------------------------------------------------
# Input

output = os.path.join(sys.argv[1], "")

# ------------------------------------------------------------

stokes_files = output+'output/phase.dat'
norm_file = output+'output/normalization.dat'
plot_dir = output+'plot/'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

phase, stokes_i, error_i, stokes_q, error_q, stokes_u, error_U, stokes_v, error_v = np.loadtxt(stokes_files, unpack=True)

if os.path.exists(norm_file):
    # Reflected light normalization (W m-2 um-1)
    _, _, norm = np.loadtxt(norm_file, unpack=True)

else:
    # Thermal light
    norm = None

stokes_pi = np.sqrt(stokes_q**2+stokes_u**2+stokes_v**2)

if norm is not None:
    # Lambertian surface
    ss_albedo = 1.
    x = np.linspace(0, 180, 1000)
    y = (2./3.)*ss_albedo*(np.sin(x*np.pi/180.)+(np.pi-(x*np.pi/180.))*np.cos(x*np.pi/180.))/np.pi
    plt.plot(x,y, ls='--')
    ss_albedo = 0.5
    x = np.linspace(0, 180, 1000)
    y = (2./3.)*ss_albedo*(np.sin(x*np.pi/180.)+(np.pi-(x*np.pi/180.))*np.cos(x*np.pi/180.))/np.pi
    plt.plot(x,y, ls='--')

# Stokes I
plt.xlabel('Phase angle (deg)')
if norm is None:
    plt.ylabel('Stokes I (W m$^{-2}$ µm$^{-1}$)')
    plt.plot(phase, stokes_i, ls='-')
else:
    plt.ylabel('Normalized Stokes I')
    plt.plot(phase, stokes_i/norm, ls='-')
    plt.ylim(0,1)
plt.xlim(0,180)
plt.savefig(plot_dir+'phase_i.pdf', bbox_inches='tight')
plt.clf()

# Stokes Q
plt.xlabel('Phase angle (deg)')
if norm is None:
    plt.ylabel('Stokes Q')
    plt.plot(phase, stokes_q, ls='-')
else:
    plt.ylabel('Normalized Stokes Q')
    plt.plot(phase, stokes_q/norm, ls='-')
    plt.ylim(-1,1)
plt.xlim(0,180)
plt.savefig(plot_dir+'phase_q.pdf', bbox_inches='tight')
plt.clf()

# Stokes U
plt.xlabel('Phase angle (deg)')
if norm is None:
    plt.ylabel('Stokes U')
    plt.plot(phase, stokes_u, ls='-')
else:
    plt.ylabel('Normalized Stokes U')
    plt.plot(phase, stokes_u/norm, ls='-')
    plt.ylim(-1,1)
plt.xlim(0,180)
plt.savefig(plot_dir+'phase_u.pdf', bbox_inches='tight')
plt.clf()

# Stokes V
plt.xlabel('Phase angle (deg)')
plt.ylabel('Normalized Stokes V')
if norm is None:
    plt.plot(phase, stokes_v, ls='-')
else:
    plt.plot(phase, stokes_v/norm, ls='-')
    plt.ylim(-1,1)
plt.xlim(0,180)
plt.savefig(plot_dir+'phase_v.pdf', bbox_inches='tight')
plt.clf()

# -Q/I
plt.xlabel('Phase angle (deg)')
plt.ylabel('-Q/I')
plt.plot(phase[stokes_i>0.], -stokes_q[stokes_i>0.]/stokes_i[stokes_i>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plot_dir+'phase_q_pol.pdf', bbox_inches='tight')
plt.clf()

# U/I
plt.xlabel('Phase angle (deg)')
plt.ylabel('U/I')
plt.plot(phase[stokes_i>0.], stokes_u[stokes_i>0.]/stokes_i[stokes_i>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plot_dir+'phase_u_pol.pdf', bbox_inches='tight')
plt.clf()

# V/I
plt.xlabel('Phase angle (deg)')
plt.ylabel('V/I')
plt.plot(phase[stokes_i>0.], stokes_v[stokes_i>0.]/stokes_i[stokes_i>0.], ls='-')
plt.xlim(0,180)
plt.ylim(-1,1)
plt.savefig(plot_dir+'phase_v_pol.pdf', bbox_inches='tight')
plt.clf()

# Degree of polarization
plt.xlabel('Phase angle (deg)')
plt.ylabel('Degree of polarization')
plt.plot(phase[stokes_i>0.], stokes_pi[stokes_i>0.]/stokes_i[stokes_i>0.], ls='-')
plt.xlim(0,180)
plt.ylim(0,1)
plt.savefig(plot_dir+'phase_polarization.pdf', bbox_inches='tight')
