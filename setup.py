from setuptools import setup

setup(
    name='seaicert',
    version='0.1.0',
    author='Andrew P. Barrett',
    author_email='andrew.barrett@colorado.edu',
    packages=["seaicert"],
    package_data={"seaicert": ["lib/libcrm.so"]},
    install_requires=[
        'pytest',
        ],
    license='LICENSE.txt',
    description='A python wrapper for the CESM2 Delta-Eddington sea ice radiative transfer code',
    long_description=open('README.md').read(),
)
