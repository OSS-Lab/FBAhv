from setuptools import setup, find_packages
from codecs     import open
from os         import path

setup(
    name            = 'FBAhv model analysis'
    version         = '1.0.0'
    description     = 'model analysis scripts for FBAhv'
    url             = 'www.osslab.lifesci.warwick.ac.uk'
    author          = 'Hadrien Dellattre & Orkun S Soyer, OSS Lab, University of Warwick'
    author_email    = 'O.Soyer@warwick.ac.uk'
    license         = 'BSD2'

    classifiers = [
        'Development Status :: 3 - Alpha'
        'Programming Language :: Python :: 3'
    ]

    keywords        = 'virus host metabolism fba cobrapy metabolic interactions infection'
    install_requires=['cobrapy'],['numpy'],['scipy']
)
