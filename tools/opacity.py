import os
import sys
import math
import shutil

import numpy as np

from astropy.io import fits
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d


def opacity_rayleigh(wavelength,
                     output="rayleigh.fits",
                     albedo=1.,
                     depolarization=0.,
                     mmw=2.02):
    """
    Function to create constant opacities with Rayleigh scattering matrix.

    :param wavelength: FITS filename that contains the wavelength points or a tuple with
                       the minimum wavelength (micron), maximum wavelength (micron), and number
                       of equally spaced wavelength points.
    :type wavelength: str or (float, float, int)
    :param output: Output FITS file name.
    :type output: str
    :param albedo: Single scattering albedo
    :type albedo: float
    :param depolarization: Depolarization factor.
    :type depolarization: float
    :param mmw: Mean molecular weight of scattering gas (g/mol).
    :type mmw: float

    :return: None
    """

    if isinstance(wavelength, str):
        hdu = fits.open(wavelength)
        wavelength = hdu[0].data[0]
        hdu.close()

    else:
        wavelength = np.linspace(wavelength[0], wavelength[1], wavelength[2], endpoint=True)

    avogadro = 6.02214129e23
    loschmidt = 2.6867805e19 # [cm-3]

    gas_mass_scat = mmw/avogadro # Molecule mass [g]

    list_wl = []
    list_scat = []
    list_abs = []
    list_ext = []

    for _, item in enumerate(wavelength):
        a = 13.58e-5
        b = 7.52e-3
        ri = 1.+a+a*b/(item**2) # Wavelength in [micron]

        rindex = (ri*ri-1.)*(ri*ri-1.)/((ri*ri+2.)*(ri*ri+2.))
        dep = (6.+3.*depolarization)/(6.-7.*depolarization)

        # Rayleigh cross section [cm2]
        cross_section = 24.*math.pi*math.pi*math.pi*rindex*dep/(((item*1.e-4)**4)*(loschmidt**2))

        scat_opacity = cross_section/gas_mass_scat

        list_wl.append(item)
        list_scat.append(scat_opacity)
        list_ext.append(scat_opacity/albedo)
        list_abs.append(scat_opacity/albedo-scat_opacity)

    opacity = np.zeros((4, len(list_wl)))
    for i in range(len(list_wl)):
        opacity[0, i] = list_wl[i]
        opacity[1, i] = list_ext[i]
        opacity[2, i] = list_abs[i]
        opacity[3, i] = list_scat[i]

    def rayleigh_p11(theta):
        alpha = math.cos(theta)
        delta = (1.-depolarization)/(1.+depolarization/2.)
        ray = (((alpha*alpha+1.)*delta)+(1.-delta))*math.sin(theta)

        return ray

    rayleigh_norm, _ = quad(rayleigh_p11, 0., math.pi)
    rayleigh_norm *= 2.*math.pi

    def rayleigh_scatter(alpha):
        matrix = np.zeros((16))

        delta = (1.-depolarization)/(1.+depolarization/2.)
        delta_prime = (1.-2.*depolarization)/(1.-depolarization)

        matrix[0] = alpha*alpha + 1.
        matrix[1] = alpha*alpha - 1.
        matrix[4] = matrix[1]
        matrix[5] = matrix[0]
        matrix[10] = 2.*alpha
        matrix[15] = delta_prime*matrix[10]

        matrix = delta * matrix
        matrix[0] = matrix[0] + (1.-delta)

        return matrix

    scatter = np.zeros((180, 16, len(list_wl)))

    for i in range(len(list_wl)):
        for j in range(180):
            matrix_low = rayleigh_scatter(math.cos(float(j)*math.pi/180.))
            matrix_up = rayleigh_scatter(math.cos(float(j+1)*math.pi/180.))

            for m in range(16):
                scatter[j, m, i] = (matrix_low[m]+matrix_up[m])/2.
                scatter[j, m, i] /= rayleigh_norm

    hdulist = fits.HDUList()
    hdulist.append(fits.ImageHDU(opacity, name='opacity'))
    hdulist.append(fits.ImageHDU(scatter, name='scatter'))
    hdu = hdulist[0]
    hdu.header['COMMENT'] = '1. Wavelength [micron]'
    hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
    hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
    hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
    hdulist.writeto(output, overwrite=True)
    hdulist.close()

def opacity_isotropic(wavelength,
                      output="isotropic.fits",
                      absorption=0.5,
                      scattering=1.0):
    """
    Function to create constant opacities with isotropic scattering matrix.

    :param wavelength: FITS filename that contains the wavelength points or a tuple with
                       the minimum wavelength (micron), maximum wavelength (micron), and number
                       of equally spaced wavelength points.
    :type wavelength: str or (float, float, int)
    :param output: Output FITS file name.
    :type output: str
    :param absorption: Absorption opacity (cm2 g-1).
    :type absorption: float
    :param scattering: Scattering opacity (cm2 g-1).
    :type scattering: float

    :return: None
    """

    if isinstance(wavelength, str):
        hdu = fits.open(wavelength)
        wavelength = hdu[0].data[0]
        hdu.close()

    else:
        wavelength = np.linspace(wavelength[0], wavelength[1], wavelength[2], endpoint=True)

    opacity = np.zeros((4, len(wavelength)))
    for i in range(len(wavelength)):
        opacity[0, i] = wavelength[i]
        opacity[1, i] = absorption + scattering
        opacity[2, i] = absorption
        opacity[3, i] = scattering

    scatter = np.zeros((180, 16, len(wavelength)))

    for i, _ in enumerate(wavelength):
        for j in range(180):
            scatter[j, 0, i] = 1./(4.*math.pi)

    hdulist = fits.HDUList()
    hdulist.append(fits.ImageHDU(opacity, name='opacity'))
    hdulist.append(fits.ImageHDU(scatter, name='scatter'))
    hdu = hdulist[0]
    hdu.header['COMMENT'] = '1. Wavelength [micron]'
    hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
    hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
    hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
    hdulist.writeto(output, overwrite=True)
    hdulist.close()

def opacity_henyey_greenstein(wavelength,
                              output="henyey.fits",
                              absorption=0.5,
                              scattering=1.0,
                              g=(0.8, 0., 0.,),
                              weight=(1., 0., 0.),
                              p_linear=0.5):
    """
    Function to create constant opacities with Henyey-Greenstein scattering matrix.

    :param wavelength: FITS filename that contains the wavelength points or a tuple with
                       the minimum wavelength (micron), maximum wavelength (micron), and number
                       of equally spaced wavelength points.
    :type wavelength: str or (float, float, int)
    :param output: Output FITS file name.
    :type output: str
    :param absorption: Absorption opacity (cm2 g-1).
    :type absorption: float
    :param scattering: Scattering opacity (cm2 g-1).
    :type scattering: float
    :param g: Asymmetry parameters of the three Henyey-Greenstein functions.
    :type g: (float, float, float)
    :param weight: Weights for the Henyey-Greenstein functions.
    :type weight: (float, float, float)
    :param p_linear: Peak of the fractional polarization.
    :type p_linear: float

    :return: None
    """

    p_circular = 0.0
    skew = 0.0

    if isinstance(wavelength, str):
        hdu = fits.open(wavelength)
        wavelength = hdu[0].data[0]
        hdu.close()

    else:
        wavelength = np.linspace(wavelength[0], wavelength[1], wavelength[2], endpoint=True)

    opacity = np.zeros((4, wavelength.shape[0]))
    for i, item in enumerate(wavelength):
        opacity[0, i] = item
        opacity[1, i] = absorption+scattering
        opacity[2, i] = absorption
        opacity[3, i] = scattering

    def hg_p11(theta):
        alpha = math.cos(theta)

        henyey = weight[0] * (1.-g[0]*g[0]) / ((1.+g[0]*g[0]-2.*g[0]*alpha)**(3./2.))
        henyey += weight[1] * (1.-g[1]*g[1]) / ((1.+g[1]*g[1]-2.*g[1]*alpha)**(3./2.))
        henyey += weight[2] * (1.-g[2]*g[2]) / ((1.+g[2]*g[2]-2.*g[2]*alpha)**(3./2.))

        return henyey*math.sin(theta)

    def hg_scatter(alpha):
        scatter_matrix = np.zeros((16))

        alpha_f = alpha * (1.+ 3.13 * skew * math.exp(-7.*alpha/math.pi))
        cos_alpha_f = math.cos(alpha_f)

        scatter_matrix[0] = weight[0] * (1.-g[0]*g[0]) / ((1.+g[0]*g[0]-2.*g[0]*alpha)**(3./2.))
        scatter_matrix[0] += weight[1] * (1.-g[1]*g[1]) / ((1.+g[1]*g[1]-2.*g[1]*alpha)**(3./2.))
        scatter_matrix[0] += weight[2] * (1.-g[2]*g[2]) / ((1.+g[2]*g[2]-2.*g[2]*alpha)**(3./2.))
        scatter_matrix[1] = -p_linear * scatter_matrix[0] * (1.-alpha*alpha) / (1.+alpha*alpha)
        scatter_matrix[4] = scatter_matrix[1]
        scatter_matrix[5] = scatter_matrix[0]
        scatter_matrix[10] = scatter_matrix[0] * (2.*alpha) / (1.+alpha*alpha)
        scatter_matrix[11] = p_circular * scatter_matrix[5] * (1.-cos_alpha_f*cos_alpha_f) / (1.+cos_alpha_f*cos_alpha_f)
        scatter_matrix[14] = -scatter_matrix[11]
        scatter_matrix[15] = scatter_matrix[10]

        return scatter_matrix

    hg_norm, _ = quad(hg_p11, 0., math.pi)
    hg_norm *= 2.*math.pi

    scatter = np.zeros((180, 16, len(wavelength)))

    for i, item in enumerate(wavelength):
        for j in range(180):
            scatter_low = hg_scatter(math.cos(float(j)*math.pi/180.))
            scatter_up = hg_scatter(math.cos(float(j+1)*math.pi/180.))

            for m in range(16):
                scatter[j, m, i] = (scatter_low[m]+scatter_up[m])/2.
                scatter[j, m, i] /= hg_norm

    hdulist = fits.HDUList()
    hdulist.append(fits.ImageHDU(opacity, name='opacity'))
    hdulist.append(fits.ImageHDU(scatter, name='scatter'))
    hdu = hdulist[0]
    hdu.header['COMMENT'] = '1. Wavelength [micron]'
    hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
    hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
    hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
    hdulist.writeto(output, overwrite=True)
    hdulist.close()

def opacity_gas(abs_file,
                wavelength=None,
                output="gas.fits",
                vmr=1.8e-3,
                mmw_abs=16.04,
                mmw_scat=2.02,
                depolarization=0.02):
    """
    Function to create wavelength dependent opacities with Rayleigh scattering matrix.

    :param abs_file: ASCII file that contains the wavelength dependent absorption coefficients.
    :type abs_file: str
    :param wavelength: Tuple with the minimum wavelength (micron), maximum wavelength (micron),
                       and number of equally spaced wavelength points. The wavelengths from the
                       absorption coefficients file are used if the number of points is set to
                       None.
    :type wavelength: (float, float, int)
    :param output: Output FITS file name.
    :type output: str
    :param vmr: Volume mixing ratio of absorbing molecule.
    :type vmr: float
    :param mmw_abs: Mean molecular weight of absorbing gas (g/mol).
    :type mmw_abs: float
    :param mmw_scat: Mean molecular weight of scattering gas (g/mol).
    :type mmw_scat: float
    :param depolarization: Depolarization factor.
    :type depolarization: float

    :return: None
    """

    if wavelength[2] is not None:
        wavelength = np.linspace(wavelength[0], wavelength[1], wavelength[2], endpoint=True)

    avogadro = 6.02214129e23 # [mol-1]
    loschmidt = 2.6867805e19 # [cm-3]

    gas_mass_abs = mmw_abs / avogadro # Molecule mass [g]
    gas_mass_scat = mmw_scat / avogadro # Molecule mass [g]

    # Wavelength [micron] - Absorption [cm2 mol-1]
    file_wl, file_abs = np.loadtxt(abs_file, unpack=True)

    file_abs = file_abs/gas_mass_abs # [cm2 molecule-1] -> [cm2 g-1]

    list_wl = []
    list_scat = []
    list_abs = []
    list_ext = []

    wavelength_new = []
    absorption_new = []

    if wavelength[2] is None:
        for i, item in enumerate(file_wl):
            if item >= wavelength[0]:
                wavelength_new.append(item)
                absorption_new.append(file_abs[i])

            if item > wavelength[1]:
                break

    else:
        y_interp = interp1d(file_wl, file_abs)

        wavelength_new = np.copy(wavelength)
        absorption_new = y_interp(wavelength_new)

    for i in range(len(wavelength_new)):
        # H2 refractive index
        a = 13.58e-5
        b = 7.52e-3 # [micron^2]
        ri = 1.+a+a*b/(wavelength_new[i]*wavelength_new[i])

        # Rayleigh cross section [cm2]
        dep = (6.+3.*depolarization)/(6.-7.*depolarization)
        rindex = ((ri*ri-1.)/loschmidt)**2.
        cross_section = (8.*math.pi**3./3.)*rindex*dep
        cross_section /= (wavelength_new[i]*1.e-4)**4.

        list_wl.append(wavelength_new[i])
        list_scat.append(cross_section/gas_mass_scat)
        list_abs.append(absorption_new[i]*vmr)
        list_ext.append(cross_section/gas_mass_scat+absorption_new[i]*vmr)

    opacity = np.zeros((4, len(list_wl)))
    for i in range(len(list_wl)):
        opacity[0, i] = list_wl[i]
        opacity[1, i] = list_ext[i]
        opacity[2, i] = list_abs[i]
        opacity[3, i] = list_scat[i]

    def rayleigh_p11(theta):
        alpha = math.cos(theta)
        delta = (1.-depolarization) / (1.+depolarization/2.)
        ray = (((alpha*alpha+1.)*delta)+(1.-delta))*math.sin(theta)

        return ray

    def rayleigh_scatter(alpha):
        delta = (1.-depolarization) / (1.+depolarization/2.)
        delta_prime = (1.-2.*depolarization) / (1.-depolarization)

        scatter_matrix = np.zeros((16))
        scatter_matrix[0] = alpha*alpha + 1.
        scatter_matrix[1] = alpha*alpha - 1.
        scatter_matrix[4] = scatter_matrix[1]
        scatter_matrix[5] = scatter_matrix[0]
        scatter_matrix[10] = 2.*alpha
        scatter_matrix[15] = delta_prime*scatter_matrix[10]
        scatter_matrix = delta * scatter_matrix
        scatter_matrix[0] = scatter_matrix[0] + (1.-delta)

        return scatter_matrix

    rayleigh_norm, _ = quad(rayleigh_p11, 0., math.pi)
    rayleigh_norm *= 2.*math.pi

    scatter = np.zeros((180, 16, len(list_wl)))
    for i in range(len(list_wl)):
        for j in range(180):
            matrix_low = rayleigh_scatter(math.cos(float(j)*math.pi/180.))
            matrix_up = rayleigh_scatter(math.cos(float(j+1)*math.pi/180.))

            for m in range(16):
                scatter[j, m, i] = (matrix_low[m]+matrix_up[m])/2.
                scatter[j, m, i] /= rayleigh_norm

    hdulist = fits.HDUList()
    hdulist.append(fits.ImageHDU(opacity, name='opacity'))
    hdulist.append(fits.ImageHDU(scatter, name='scatter'))
    hdu = hdulist[0]
    hdu.header['COMMENT'] = '1. Wavelength [micron]'
    hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
    hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
    hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
    hdulist.writeto(output, overwrite=True)
    hdulist.close()

def opacity_dhs(ri_file,
                wavelength=None,
                output="dhs.fits",
                nr=100,
                nf=20,
                density=1.,
                amin=0.1,
                amax=1.0,
                apow=3.5,
                fmax=0.,
                r_eff=1.0,
                v_eff=0.1):
    """
    Function to create DHS opacities and scattering matrices.

    :param ri_file: ASCII file that contains the wavelength dependent complex refractive index.
    :type ri_file: str
    :param wavelength: FITS filename that contains the wavelength points or a tuple with
                       the minimum wavelength (micron), maximum wavelength (micron), and number
                       of equally spaced wavelength points. The wavelengths from the refractive
                       index file are used if the number of points is set to None.
    :type wavelength: str or (float, float, int)
    :param output: Output FITS file name.
    :type output: str
    :param nr: Number of radius points.
    :type nr: int
    :param nf: Number of volume fractions for distribution of hollow spheres (DHS).
    :type nf: int
    :param density: Particle density (g cm-3).
    :type density: float
    :param amin: Minimum particle size (micron).
    :type amin: float
    :param amax: Maximum particle size (micron).
    :type amax: float
    :param apow: Size distribution power law index
    :type apow: float
    :param fmax: Irregularity parameter, maximum volume void fraction for DHS. Setting fmax to 0
                 equals Mie theory.
    :type fmax: float
    :param r_eff: Effective radius (micron), overrules amin, amax, apow.
    :type r_eff: float
    :param v_eff: Effective variance (dimensionless), overrules amin, amax, apow.
    :type v_eff: float

    :return: None
    """

    code_path = os.path.dirname(os.path.abspath(__file__))

    # Percentage set to 100% (i.e., no mixture of multiple grain species)
    percentage = 100.

    if sys.platform == "linux" or sys.platform == "linux2":
        computepart = code_path[:-6]+'/bin/computepart_linux'
    elif sys.platform == "darwin":
        computepart = code_path[:-6]+'/bin/computepart_mac'

    temp_dir = "temp/"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    os.chdir(temp_dir)

    file_wl, _, _ = np.loadtxt(ri_file, unpack=True)

    f = open("wavelength.dat", 'w')

    if isinstance(wavelength, str):
        hdu = fits.open("../"+wavelength)
        wavelength = hdu[0].data[0]
        hdu.close()

        for i in range(len(wavelength)):
            f.write(str(wavelength[i])+"\n")

    elif wavelength[2] is None:
        if np.size(file_wl) == 1:
            f.write(str(file_wl))

        else:
            for i, item in enumerate(file_wl):
                if item >= wavelength[0] and item <= wavelength[1]:
                    f.write(str(file_wl[i])+"\n")

    else:
        wavelength = np.linspace(wavelength[0], wavelength[1], wavelength[2], endpoint=True)

        for _, item in enumerate(wavelength):
            f.write(str(item)+"\n")

    f.close()

    f = open("mie.in", 'w')
    f.write(str(nr)+"\n")
    f.write(str(nf)+"\n")
    f.write("'"+str(ri_file)+"' \n")
    f.write(str(percentage)+"\t"+str(density)+"\t"+str(amin)+"\t"+str(amax)+"\t"+str(apow)+"\t"+str(fmax))
    f.close()

    os.chmod(computepart, 700)

    if r_eff > 0.:
        os.system(computepart+" mie.in wavelength.dat "+str(r_eff)+" "+str(v_eff))
    else:
        os.system(computepart+" mie.in wavelength.dat")

    os.chdir("../")
    os.rename(temp_dir+"particle.fits", output)
    shutil.rmtree(temp_dir)

    hdulist = fits.open(output)
    wavelengths = int(hdulist[1].header['NAXIS1'])
    elements = int(hdulist[1].header['NAXIS2'])
    opacity = hdulist[0].data
    scatter = hdulist[1].data
    hdulist.close()

    # Change to 16-elements scatter matrix

    scatterNew = np.zeros((180, 16, wavelengths))

    for i in range(wavelengths):
        for j in range(180):
            scatterNew[j, 0, i] = scatter[j, 0, i]
            scatterNew[j, 1, i] = scatter[j, 1, i]
            scatterNew[j, 4, i] = scatter[j, 1, i]
            scatterNew[j, 5, i] = scatter[j, 2, i]
            scatterNew[j, 10, i] = scatter[j, 3, i]
            scatterNew[j, 11, i] = scatter[j, 4, i]
            scatterNew[j, 14, i] = -scatter[j, 4, i]
            scatterNew[j, 15, i] = scatter[j, 5, i]

    scatter = scatterNew

    # Scatter matrix normalization

    angle = np.zeros(180)
    for i in range(180):
        angle[i] = (float(i)+0.5) * math.pi/180.

    for j in range(wavelengths):

        p11Int = scatter[:, 0, j]
        norm = simps(p11Int*np.sin(angle), angle)
        norm *= 2.*math.pi
        scatter[:, :, j] /= norm

    # Write updated FITS file

    hdulist = fits.HDUList()
    hdulist.append(fits.ImageHDU(opacity, name='opacity'))
    hdulist.append(fits.ImageHDU(scatter, name='scatter'))
    hdu = hdulist[0]
    hdu.header['COMMENT'] = '1. Wavelength [micron]'
    hdu.header['COMMENT'] = '2. Extinction [cm2 g-1]'
    hdu.header['COMMENT'] = '3. Absorption [cm2 g-1]'
    hdu.header['COMMENT'] = '4. Scattering [cm2 g-1]'
    hdulist.writeto(output, overwrite=True)
    hdulist.close()

def opacity_molecules(pt_file,
                      wavelength=(0.5, 20.0),
                      mmw=2.3,
                      depolarization=0.0):
    """
    Function to create opacities and scattering matrices for molecules.

    :param pt_file: ASCII file with the pressure-temperature profile.
    :type pt_file: str
    :param wavelength: Tuple with the minimum and maximum wavelength (micron)
    :type wavelength: (float, float)
    :param depolarization: Depolarization factor.
    :type depolarization: float
    :param mmw: Mean molecular weight (g/mol).
    :type mmw: float

    :return: None
    """

    dat_path = os.path.dirname(os.path.abspath(__file__))[:-6]+'/dat/molecules/'

    if not os.path.exists('opacity/'):
        os.makedirs('opacity/')

    pressure, temperature = np.loadtxt(pt_file, unpack=True)

    def get_opacity(filenumber):
        """
        :param filenumber: Filenumber corresponds with number in PTgrid.dat.
        :type filenumber: int

        :return: Wavelength and opacity array.
        :rtype: ndarray, ndarray
        """

        wl, op = np.loadtxt(dat_path+'opacity_aver_'+str(int(filenumber)).zfill(4)+'.dat', unpack=True)

        return wl, op

    def get_pt(log_pressure_layer,
               temp_layer):
        """
        :param log_pressure_layer: Logarithm of the pressure (bar).
        :type log_pressure_layer: float
        :param temp_layer: Temperature (K).
        :type temp_layer: float

        :return: Indices array with the pressure and temperature boundaries in PTgrid.dat.
                 Note that indices start at zero and filenames at 1.
        :rtype: ndarray
        """

        opacity_pt_list = np.genfromtxt(dat_path+'PTgrid.dat', skip_header=1)
        pressure_layer = 10.**log_pressure_layer

        index_opac = opacity_pt_list[::, 0]
        p_opac = opacity_pt_list[::, 1]
        t_opac = opacity_pt_list[::, 2]

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

            return [upperPupperTindex, lowerPupperTindex, upperPlowerTindex, lowerPlowerTindex]

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

            x = p_opac[indices[0]], p_opac[indices[1]], p_opac[indices[2]], p_opac[indices[3]]
            y = t_opac[indices[0]], t_opac[indices[1]], t_opac[indices[2]], t_opac[indices[3]]

            return indices

        except UnboundLocalError:
            return [0, 0, 0, 0]

    def interpolate_opacity(pressure_layer,
                            temp_layer,
                            indices,
                            opacity_array):
        """
        :param pressure_layer: Pressure (bar).
        :type pressure_layer: float
        :param temp_layer: Temperature (K).
        :type temp_layer: float
        :param indices: Indices array of P/T boundaries in PTgrid.dat.
        :type indices: ndarray
        :param opacity_array: Array of opacities for each PT point in the 4 indexed points.
        :type opacity_array: ndarray

        :return: Opacities.
        :rtype: ndarray
        """

        opacity_pt_list = np.genfromtxt(dat_path+'PTgrid.dat', skip_header=1)

        p1 = opacity_pt_list[indices[1]][1]
        p2 = opacity_pt_list[indices[0]][1]
        t1 = opacity_pt_list[indices[2]][2]
        t2 = opacity_pt_list[indices[0]][2]

        p1 = np.log10(p1)
        p2 = np.log10(p2)
        t1 = np.log10(t1)
        t2 = np.log10(t2)

        opacity_array = np.log10(opacity_array)
        temp_layer = np.log10(temp_layer)
        pressure_layer = np.log10(pressure_layer)

        opacity_array[opacity_array < -500] = -500

        if ((p1 == p2) and (t1 == t2)):
            return 10.**opacity_array[0]

        elif (p1 == p2):
            opacities_interp = opacity_array[2] + (opacity_array[0] -opacity_array[2]) * (temp_layer-t1)/(t2-t1)

            return 10.**opacities_interp

        elif (t1 == t2):
            opacities_interp =  opacity_array[1] + (opacity_array[0] -opacity_array[1]) * (pressure_layer-p1)/(p2-p1)

            return 10.**opacities_interp

        rr1 = ((p2-pressure_layer) / (p2-p1)) * opacity_array[3] + ((pressure_layer-p1)/(p2-p1)) * opacity_array[2]
        rr2 = ((p2-pressure_layer) / (p2-p1)) * opacity_array[1] + ((pressure_layer-p1)/(p2-p1)) * opacity_array[0]
        opacities_interp = ((t2-temp_layer) / (t2-t1)) * rr1 + ((temp_layer-t1)/(t2-t1)) * rr2

        return 10.**opacities_interp

    # Opacities

    pressure_log = np.log10(pressure)

    for i in range(len(pressure_log)):
        indices = get_pt(pressure_log[i], temperature[i])

        # Just to count the number of wavelengths
        wavelength, opacity = get_opacity(indices[0]+1)

        opacityNew = np.zeros(len(wavelength))

        wavelength, opacity1 = get_opacity(indices[0]+1)
        wavelength, opacity2 = get_opacity(indices[1]+1)
        wavelength, opacity3 = get_opacity(indices[2]+1)
        wavelength, opacity4 = get_opacity(indices[3]+1)

        opacity = [opacity1, opacity2, opacity3, opacity4]
        opacityNew = interpolate_opacity(pressure[i], temperature[i], indices, opacity)

        # Individual molecule opacity and total absorption for each layer
        absFile = "opacity/absorption_"+str('%02d' % (len(pressure_log)-i))+".dat"
        f = open(absFile, 'w')
        f.write("# Wavelength [micron] - Opacity x vmr [cm2/molecule]\n\n")
        for i in range(len(wavelength)):
            f.write(str(wavelength[i])+'\t'+str(opacityNew[i])+'\n')
        f.close()

    # Make FITS files

    opacityDir = 'opacity/'

    avogadro = 6.02214129e23
    loschmidt = 2.6867805e19 # [cm-3]

    def refractiveIndexH2(wavelength):
        a = 13.58e-5
        b = 7.52e-3
        ri = 1. + a + a*b / (wavelength*wavelength)

        return ri

    def rayleigh_p11(theta):
        alpha = math.cos(theta)
        delta = (1.-depolarization) / (1.+depolarization/2.)
        ray = (((alpha*alpha+1.)*delta)+(1.-delta))*math.sin(theta)

        return ray

    def rayleigh_scatter(alpha):
        delta = (1.-depolarization) / (1.+depolarization/2.)
        delta_prime = (1.-2.*depolarization) / (1.-depolarization)

        scatter_matrix = np.zeros((16))
        scatter_matrix[0] = alpha*alpha + 1.
        scatter_matrix[1] = alpha*alpha - 1.
        scatter_matrix[4] = scatter_matrix[1]
        scatter_matrix[5] = scatter_matrix[0]
        scatter_matrix[10] = 2.*alpha
        scatter_matrix[15] = delta_prime*scatter_matrix[10]
        scatter_matrix = delta * scatter_matrix
        scatter_matrix[0] = scatter_matrix[0] + (1.-delta)

        return scatter_matrix

    rayleigh_norm, _ = quad(rayleigh_p11, 0., math.pi)
    rayleigh_norm *= 2.*math.pi

    for file in os.listdir(opacityDir):
        if file.startswith('absorption_') and file.endswith('.dat'):
            wavelength, absorption = np.loadtxt(opacityDir+file, unpack=True)
            mass = mmw / avogadro # Molecule mass [g]
            absorption /= mass # [cm2 molecule-1] -> [cm2 g-1]

            scatter = np.zeros((180, 16, len(wavelength)))

            list_wl = []
            list_scat = []
            list_abs = []
            list_ext = []

            # def rayleigh_h2(wavelength):
            #     sigma = (8.14e-13/(wavelength**4)) + (1.28e-6/(wavelength**6)) + (1.61/(wavelength**8))
            #     return sigma

            for i in range(len(wavelength)):
                if wavelength[i] >= wavelength[0]:
                    # Rayleigh cross section [cm2]
                    ri = refractiveIndexH2(wavelength[i])
                    rindex = (ri*ri-1.)*(ri*ri-1.)/((ri*ri+2.)*(ri*ri+2.))
                    dep = (6.+3.*depolarization)/(6.-7.*depolarization)
                    crossSection = 24.*math.pi*math.pi*math.pi*rindex*dep/(((wavelength[i]*1.e-4)**4)*(loschmidt**2))

                    list_wl.append(wavelength[i])
                    list_scat.append(crossSection/mass)
                    list_abs.append(absorption[i])
                    list_ext.append(crossSection/mass+absorption[i])

                    if wavelength[i] > wavelength[1]:
                        break

            opacity = np.zeros((4, len(list_wl)))
            scatter = np.zeros((180, 16, len(list_wl)))

            for i in range(len(list_wl)):
                opacity[0, i] = list_wl[i] # [micron]
                opacity[1, i] = list_ext[i] # [cm2 g-1]
                opacity[2, i] = list_abs[i] # [cm2 g-1]
                opacity[3, i] = list_scat[i] # [cm2 g-1]

            # Scattering matrices

            scatter = np.zeros((180, 16, len(list_wl)))

            for i in range(len(list_wl)):
                for j in range(180):
                    scatter_low = rayleigh_scatter(math.cos(float(j)*math.pi/180.))
                    scatter_up = rayleigh_scatter(math.cos(float(j+1)*math.pi/180.))

                    for m in range(16):
                        scatter[j, m, i] = (scatter_low[m]+scatter_up[m])/2.
                        scatter[j, m, i] /= rayleigh_norm

            fitsfile = opacityDir+'gas_opacity_'+file[11:13]+'.fits'

            hdulist = fits.HDUList()
            hdulist.append(fits.ImageHDU(opacity, name='opacity'))
            hdulist.append(fits.ImageHDU(scatter, name='scatter'))
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
