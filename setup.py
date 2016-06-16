#!/usr/bin/evn python

import sys
import os
from setuptools import setup, find_packages

setup(name='gtapps_mp',
      version='1.5',
      description='Fermi LAT Multicore Scripts',
      author='Jeremy S. Perkins (FSSC)',
      author_email='fermihelp@milkyway.gsfc.nasa.gov',
      url='http://fermi.gsfc.nasa.gov/ssc/',
      packages=find_packages(exclude='tests'),
      entry_points= {'console_scripts': [
                      'gtdiffrsp_mp = gtapps_mp.gtdiffrsp_mp:cli',
                      'gtexpmap_mp = gtapps_mp.gtexpmap_mp:cli',
                      'gtltcube_mp = gtapps_mp.gtltcube_mp:cli',
                      'gttsmap_mp = gtapps_mp.gttsmap_mp:cli',
                      ]},
      )
