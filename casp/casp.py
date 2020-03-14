import csv

from scipy import interpolate
from scipy.integrate import quad

from casp.utils import *


def prob_bloom(R_frb, m_i, R_0=0.2, R_h=0.25):
    """

    :param R_frb: The 1 sigma localization radius of the FRB
    :param m_i: r band magnitude of the galaxy
    :param R_0: Radial angular separation between the FRB position and a presumed host
    :param R_h: Galaxy half light radius
    :return: Probability of chance coincidence
    """
    r_i = get_R(R_frb, R_0, R_h)
    factor = 3600 ** 2 * 0.334 * np.log(10)
    # galaxy per arcsecond square
    mean_surfden_gal = (1 / factor) * 10 ** (0.334 * (m_i - 22.963) + 4.320)
    num_gal = np.pi * (r_i ** 2) * mean_surfden_gal
    return 1 - np.exp(-1 * num_gal)


def prob_eb17(R_frb, m, R_0=0.2, R_h=0.25, ret_numgal=False):
    """

    :param R_frb: The 1 sigma localization radius of the FRB
    :param m: r band magnitude of the galaxy
    :param R_0: Radial angular separation between the FRB position and a presumed host
    :param R_h: Galaxy half light radius
    :param ret_numgal: to return the number of galaxies along with the chance coincidence probability
    :return: Probability of chance coincidence
    """
    r_dat, mag_uniq, cvs = read_r_mags()

    spl = interpolate.UnivariateSpline(x=mag_uniq,
                                       y=np.log10(r_dat),
                                       bbox=[-100, 100],
                                       k=3)

    def n_gal(m_r):
        return 10 ** spl(m_r)

    num_dens_gal = quad(n_gal, 0, m)[0]
    R = get_R(R_frb, R_0, R_h)

    deg2arcsec = 60 * 60
    num_gals = np.pi * (R / deg2arcsec) ** 2 * num_dens_gal

    if ret_numgal:
        return 1 - np.exp(-1 * num_gals), num_gals
    else:
        return 1 - np.exp(-1 * num_gals)


def prob_eb17_z(z, R_frb, R_0=0.2, R_h=0.25, ret_numgal=False, num_galaxies_file=None):
    """

    :param z: Maximum redshift of the FRB
    :param R_frb: The 1 sigma localization radius of the FRB
    :param R_0: Radial angular separation between the FRB position and a presumed host
    :param R_h: Galaxy half light radius
    :param ret_numgal: To return number of galaxies
    :param gal_num_file: File with number of galaxies in diff. redshift bins
    :return: Probability of chance coincidence
    """
    if num_galaxies_file:
        with open(num_galaxies_file, newline='') as file:
            reader = csv.reader(file, quoting=csv.QUOTE_NONNUMERIC)
            z_mins = []
            z_maxs = []
            num_total = []
            for row in reader:
                z_mins.append(row[0])
                z_maxs.append(row[1])
                num_total.append(row[2])
            z_mins = np.array(z_mins)
            z_maxs = np.array(z_maxs)
    else:
        z_mins, z_maxs, num_total = calc_num_of_galaxies()

    if z <= 1.2:
        idx = np.where((z >= z_mins) & (z < z_maxs))[0][0]
        z_max = z_maxs[idx]
        z_min = z_mins[idx]
    else:
        print('Cannot calculate for z > 1.2')
        return -1

    n = 0
    for i in range(idx):
        n += num_total[i]

    R = get_R(R_frb, R_0, R_h)
    f_a = np.pi * R ** 2 / (5.346 * 10 ** 11)
    num_gals = f_a * n

    if ret_numgal:
        return 1 - np.exp(-1 * num_gals), num_gals
    else:
        return 1 - np.exp(-1 * num_gals)
