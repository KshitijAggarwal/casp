# CASP

[![DOI](https://zenodo.org/badge/246733557.svg)](https://zenodo.org/badge/latestdoi/246733557)


Calculating Asociation Probability of FRBs

If you just want a quick estimate for a few sources, you can use the web tool [here](https://kshitijaggarwal.github.io/casp/).

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
    Chance coincidence probability (Eftekhari et al) using redshift: 0.0095

Also, check out the example notebook ([here](https://github.com/KshitijAggarwal/casp/blob/master/examples/eb17_plots.ipynb)), where I have tried to reproduce some figures from [Eftekhari et al 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...849..162E/abstract).

Requirements
---
* astropy
* numpy
* scipy
 
Citation
---
Please cite the following papers if you used `casp`: 
* Aggarwal et al (in prep.)
* [Eftekhari et al 2017](https://ui.adsabs.harvard.edu/abs/2017ApJ...849..162E/abstract)
* [Bloom et al 2003](https://ui.adsabs.harvard.edu/abs/2002AJ....123.1111B/abstract)

