#!/usr/bin/evn python

import sys
import os
from setuptools import setup, find_packages

fermi_dir = os.environ.get("FERMI_DIR")

if fermi_dir:
    print "The Fermi Science tools seem to be set up."
    print "Will install into "+fermi_dir+"."
else:
    print "The Fermi Science tools are not set up."
    sys.exit()

setup(name='gtapps_mp',
      version='1.3',
      description='Fermi LAT Multicore Scripts',
      author='Jeremy S. Perkins (FSSC)',
      author_email='fermihelp@milkyway.gsfc.nasa.gov',
      url='http://fermi.gsfc.nasa.gov/ssc/',
      packages=find_packages(exclude='tests'),
      scripts=['scripts/gtdiffrsp_mp',
               'scripts/gtexpmap_mp',
               'scripts/gtltcube_mp',
               'scripts/gttsmap_mp'],
#      py_modules=['gtdiffrsp_mp',
#                  'gtexpmap_mp',
#                  'gtltcube_mp',
#                  'gttsmap_mp'],
#      data_files=[(fermi_dir+"/bin",['scripts/gtdiffrsp_mp',
#                                     'scripts/gtexpmap_mp',
#                                     'scripts/gtltcube_mp',
#                                     'scripts/gttsmap_mp'])],
      )
