# CASP
Calculating Asociation Probability of FRBs

Install
---
    git clone https://github.com/KshitijAggarwal/casp.git
    cd casp
    python setup.py install

The installation will put `calc_prob.py` in your `PYTHONPATH`.

Usage
---
    calc_prob.py -r 1 -m 22 -z 0.2

Output will look like:

    Chance coincidence probability (Bloom et al): 0.0125
    Chance coincidence probability (Eftekhari et al): 0.0036
    Chance coincidence probability (Eftakhari et al) using redshift: 0.0015

Requirements
---
* astropy
* numpy
* scipy
 
