#!/usr/bin/env python3

import argparse, os
from casp import casp

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate association probability of FRBs.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--frb_loc_radius', type=float, help='Localisation radius (arcsec) of FRB',
                        required=True)
    parser.add_argument('-ro', '--r_o', type=float, help='Radial angular separation between the FRB position '
                                                      'and a presumed host', required=False, default=0.2)
    parser.add_argument('-rh', '--r_h', type=float, help='Galaxy half light radius', required=False, default=0.25)
    parser.add_argument('-m', '--mag', type=float, help='r-band magnitude of the source', required=False)
    parser.add_argument('-z', '--max_redshift', help='maximum redshift of the FRB', type=float, required=False)
    values = parser.parse_args()

    if values.mag:
        p_bloom = casp.prob_bloom(R_frb=values.frb_loc_radius, m_i=values.mag, R_0=values.r_o, R_h=values.r_h)
        p_eb17 = casp.prob_eb17(R_frb=values.frb_loc_radius, m=values.mag, R_0=values.r_o, R_h=values.r_h, ret_numgal=False)
        print(f'Chance coincidence probability (Bloom et al): {p_bloom:.4f}')
        print(f'Chance coincidence probability (Eftekhari et al): {p_eb17:.4f}')

    if values.max_redshift:
        file = '/'.join(__file__.split('/')[:-3])+'/casp/data/num_galaxies.txt'
        if not os.path.isfile(file):
            file = None
        p_eb17_z = casp.prob_eb17_z(z=values.max_redshift, R_frb=values.frb_loc_radius, R_0=values.r_o, R_h=values.r_h,
                                    ret_numgal=False, num_galaxies_file=file)
        print(f'Chance coincidence probability (Eftekhari et al) using redshift: {p_eb17_z:.4f}')
