#!/usr/bin/env python3

import argparse
from casp import casp

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate association probability of FRBs.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--frb_loc_radius', type=float, help='Localisation radius (arcsec) of FRB',
                        required=True)
    parser.add_argument('-ro', '--r_o', type=float, help='Radial angular separation between the FRB position '
                                                      'and a presumed host', required=False, default=0.2)
    parser.add_argument('-rh', '--r_h', type=float, help='Galaxy half light radius', required=False, default=0.25)
    parser.add_argument('-m', '--mag', type=float, help='r-band magnitude of the source', required=True)
    parser.add_argument('-z', '--max_redshift', help='maximum redshift of the FRB', type=float, required=False)
    values = parser.parse_args()

    p_bloom = casp.prob_bloom(R_frb=values.r, m_i=values.m, R_0=values.ro, R_h=values.rh)
    p_eb17 = casp.prob_eb17(R_frb=values.r, m=values.m, R_0=values.ro, R_h=values.rh, ret_numgal=False)

    print(f'Chance coincidence probability (Bloom et al): {p_bloom}')
    print(f'Chance coincidence probability (Eftekhari et al): {p_eb17}')

    if values.z:
        p_eb17_z = casp.prob_eb17_z(z=values.z, R_frb=values.r, R_0=values.ro, R_h=values.rh, ret_numgal=False,
                    num_galaxies_file=None)
        print(f'Chance coincidence probability (Eftakhari et al) using redshift: {p_eb17_z}')
