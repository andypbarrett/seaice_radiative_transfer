from setuptools import setup

setup(
    name='seaice_radiative_transfer',
    version='0.1.0',
    author='Andrew P. Barrett',
    author_email='andrew.barrett@colorado.edu',
    packages=["seaicert"],
    package_data={"seaicert": ["1D_dE_CCSM/libcrm.so"]},
    license='LICENSE.txt',
    description='A python wrapper for the CESM2 Delta-Eddington sea ice radiative transfer code',
    long_description=open('README.md').read(),
)
