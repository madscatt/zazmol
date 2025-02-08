'''
    SASMOL  Copyright (C) 2009-2016 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os
##import setuptools
#from numpy.distutils.core import Extension, setup
from setuptools import Extension, setup
#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/05/2015      --      for sasmol distribution     :       jc
#       08/08/2016      --      branched for refactoring    :       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      setup.py is the script to install and/or update sasmol

	> sudo python setup.py install

'''

import os
import numpy

os.environ['NPY_DISABLE_CPU_FEATURES'] = 'AVX512FP16'

NUMPY_INCLUDE = numpy.get_include()

def read(fname):
    ''' method to read a file name '''
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='sasmol',
      version='2.0',
      author='Joseph E. Curtis',
      author_email='joseph.curtis@nist.gov',
      license='GPL 3',
      url='https://github.com/madscatt/sasmol',
      platforms='Linux, Mac OS X',
      description=("A library of methods to write software for molecular systems"),
      long_description=read('README.md'),
      classifiers=["Development Status :: 4 - Beta",
                   "License :: OSI Approved :: GNU Public License 3",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Operating System :: Linux :: MacOS :: MacOS X",
                   "Programming Language :: Python :: C :: Fortran",
                   "Topic :: Scientific/Engineering :: Chemistry :: Physics"],

      package_dir={'sasmol':os.path.join('src','python')},

      packages=['sasmol','sasmol.test_sasmol','sasmol.test_sasmol.utilities','sasmol.test_sasmol.manual_tests','sasmol.test_sasmol.data','sasmol.test_sasmol.data.pdb_common','sasmol.test_sasmol.data.dcd_common','sasmol.test_sasmol.data.sasmol','sasmol.test_sasmol.data.sasmol.calculate','sasmol.test_sasmol.data.sasmol.file_io','sasmol.test_sasmol.data.sasmol.file_io.test-results','sasmol.test_sasmol.data.sasmol.linear_algebra','sasmol.test_sasmol.data.sasmol.system','sasmol.test_sasmol.data.sasmol.operate','sasmol.test_sasmol.data.sasmol.properties','sasmol.test_sasmol.test_calculate','sasmol.test_sasmol.test_file_io','sasmol.test_sasmol.test_linear_algebra','sasmol.test_sasmol.test_operate','sasmol.test_sasmol.test_properties','sasmol.test_sasmol.test_subset','sasmol.extensions','sasmol.extensions.dcdio','sasmol.extensions.mask','sasmol.extensions.matrix_math'],
#      packages=['sasmol','sasmol.test_sasmol','sasmol.test_sasmol.utilities','sasmol.test_sasmol.manual_tests','sasmol.test_sasmol.data','sasmol.test_sasmol.data.pdb_common','sasmol.test_sasmol.data.dcd_common','sasmol.test_sasmol.data.sasmol','sasmol.test_sasmol.data.sasmol.calculate','sasmol.test_sasmol.data.sasmol.file_io','sasmol.test_sasmol.data.sasmol.file_io.test-results','sasmol.test_sasmol.data.sasmol.linear_algebra','sasmol.test_sasmol.data.sasmol.system','sasmol.test_sasmol.data.sasmol.operate','sasmol.test_sasmol.data.sasmol.properties','sasmol.test_sasmol.test_calculate','sasmol.test_sasmol.test_file_io','sasmol.test_sasmol.test_linear_algebra','sasmol.test_sasmol.test_operate','sasmol.test_sasmol.test_properties','sasmol.test_sasmol.test_subset','sasmol.extensions','sasmol.extensions.dcdio','sasmol.extensions.view','sasmol.extensions.mask','sasmol.extensions.matrix_math'],

    ext_modules=[
    Extension('sasmol._dcdio',
                   sources=[os.path.join('src','python','extensions','dcdio','dcdio.c'),
                            os.path.join('src','python','extensions','dcdio','dcdio_module.c')],
                   include_dirs=[numpy.get_include()]),

#    Extension('sasmol._view_vmd',[os.path.join('src','python','extensions','view','view_vmd.i'),os.path.join('src','python','extensions','view','view_vmd.c'),os.path.join('src','python','extensions','view','imd.c'),os.path.join('src','python','extensions','view','vmdsock.c')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol._mask',[os.path.join('src','python','extensions','mask','mask.i'),os.path.join('src','python','extensions','mask','mask.c')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol.overlap',[os.path.join('src','python','extensions','overlap','overlap.c')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol.matrix_math',[os.path.join('src','python','extensions','matrix_math','matrix_math.c')],include_dirs=[NUMPY_INCLUDE])
    ],

      data_files=[
          (os.path.join('src','python'), [os.path.join('src', 'python', 'extensions', 'mask', 'mask.py')])
      ]

    )

