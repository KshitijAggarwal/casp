from setuptools import setup

setup(
    name='casp',
    version='0.1.0',
    packages=['casp'],
    scripts=['bin/calc_prob.py'],
    url='https://github.com/KshitijAggarwal/casp',
    author='Kshitij Aggarwal',
    author_email='ka0064@mix.wvu.edu',
    license='',
    description='Calculating Association Probability of FRBs',
    install_requires=['astropy', 'numpy', 'scipy'],
    include_package_data=True,
    package_data={
        "casp": ["data/*.txt",
                 "data/Table3MRT.fits"],
    }
)
