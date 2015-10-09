#! /usr/bin/env/python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import matplotlib
import os
import shutil


config = {
    'description': 'CMIP5 data analysis on JASMIN',
    'author': 'Angus Ferraro',
    'url': 'None',
    'author_email': 'a.j.ferraro@exeter.ac.uk.',
    'version': '0.1',
    'install_requires': ['iris'],
    'extras_require':{'testing':['nose']},
    'packages': ['jazz'],
    'name': 'jazz'
}

setup(**config)

# Install matplotlib styles

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

mpl_dir = '{0}/stylelib/'.format(matplotlib.get_configdir())
ensure_dir(mpl_dir)
shutil.copy('./data/mpl_styles/jazz_paper.mplstyle', mpl_dir)
shutil.copy('./data/mpl_styles/jazz_poster.mplstyle', mpl_dir)
shutil.copy('./data/mpl_styles/jazz_presentation.mplstyle', mpl_dir)

