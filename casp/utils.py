from astropy.cosmology import WMAP9 as cosmo
from scipy.integrate import quad
from astropy.io import fits
import numpy as np
import os


def get_R(R_frb, R_0=0.2, R_h=0.25):
    """
    Calculates Radius of localisation region in arcsecond
    Based on Bloom et al 2002 and Eftakhari et al 2017

    :param R_frb: The 1 sigma localization radius of the FRB
    :param R_0: Radial angular separation between the FRB position and a presumed host
    :param R_h: Galaxy half light radius
    :return: radius
    """
    return np.max([2*R_frb, np.sqrt(R_0**2 + 4*R_h**2)])


def read_r_mags():
    """
    Reading data used in Driver et al (2016).
    Returns R band magnitudes, Magnitude bins
    and Cosmic Variance

    :return:
    """
    data_path = os.path.split(__file__)[0] + '/data'
    table = fits.open(data_path+'/Table3MRT.fits')
    data = table[1].data
    r_mask = np.concatenate((np.where(data['Filtername'] == 'r')[0],
                             np.where(data['Filtername'] == 'F606W')[0]))
    _magbin = data['MagBinCentre'][r_mask]
    args = _magbin.argsort()
    magbin = _magbin[args]

    r_band_data = data['N(m)'][r_mask][args]
    cv = data['CosmicVariance'][r_mask][args]
    mag_uniq = np.unique(magbin)

    cvs = []
    r_dat = []
    for mag in mag_uniq:
        loc = np.where(magbin == mag)[0]
        if len(loc) > 1:
            cv_s = cv[loc]
            min_cv_loc = np.where(cv_s == cv_s.min())
            r_dat.append(r_band_data[loc[min_cv_loc]][0])
            cvs.append(cv_s[min_cv_loc][0])
        else:
            cvs.append(cv[loc][0])
            r_dat.append(r_band_data[loc][0])
    cvs = np.array(cvs)
    r_dat = np.array(r_dat)
    
    return r_dat, mag_uniq, cvs


def schechter_fn(M, phi_star, M_star, alpha):
    """
    Schechter Function

    :param M: Absolute magnitude of a given waveband (M - 5*logh)
    :param phi_star: Normalizing factor (units: 10^(-3)*h^3*Mpc^(-3)*mag^(-1))
    :param M_star: Transition from a power law luminosity function to exponential (M - 5*logh)
    :param alpha: Slope of power law variation at faint end
    :return: phi_m: (units: 10^(-3)*h^3*Mpc^(-3)*mag^(-1))
    """
    factor = -0.4*(M - M_star)
    phi_m = 0.4*np.log(10)*phi_star*(10**(factor*(alpha + 1)))*(np.exp(-10**factor))
    return phi_m


def calc_num_of_galaxies(M_min=-24, M_max=None, save = False):
    """
    :param M_min: Minimum absolute magnitude
    :param M_max: Maximum absolute magnitude
    :param save: To save the galaxy numbers to a file
    :return: redshift bin edges and number of galaxies
    """
    alphas = [-1.05, -1.175, -1.3, -1.3, -1.3, -1.3, -1.3]
    phi_stars = [14.9*0.7**3, 4.32, 3.53, 3.23, 3.46, 3.78, 2.5] # 10^-3 h^3 Mpc^-3 mag^-1  (h = 1)
    Mb_stars = [-20.44, -20.54, -20.64, -20.97, -21.08, -21.22, -21.39]  # + 5 logh70         
    
    z_mins = np.array([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
    z_maxs = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    
    num_total = []
    for idx, (z_min, z_max) in enumerate(zip(z_mins, z_maxs)):
        phi_star = phi_stars[idx]
        Mb_star = Mb_stars[idx]
        if not M_max:
            M_max = Mb_star + 5
        alpha = alphas[idx]

        V_c_max = cosmo.comoving_volume(z_max).value # Mpc^3
        V_c_min = cosmo.comoving_volume(z_min).value # Mpc^3
        num_dens_gal = quad(schechter_fn, M_min, M_max, args=(phi_star, Mb_star, alpha))[0]*10**(-3) # Mpc^-3

        num_total.append(num_dens_gal*(V_c_max - V_c_min))
        
    if save:
        data_path = os.path.split(__file__)[0] + '/data'
        with open(f"{data_path}/num_galaxies.txt","w") as f:
            for (z_min, z_max, num_gal) in zip(z_mins, z_maxs, num_total):
                f.write("{0},{1},{2}\n".format(z_min, z_max, num_gal))
                
    return z_mins, z_maxs, num_total