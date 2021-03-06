'''
    SASMOL  Copyright (C) 2009-2016 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os
#import setuptools
from numpy.distutils.core import Extension, setup
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

import numpy

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

      packages=['sasmol','sasmol.test_sasmol','sasmol.test_sasmol.utilities','sasmol.test_sasmol.manual_tests','sasmol.test_sasmol.data','sasmol.test_sasmol.data.pdb_common','sasmol.test_sasmol.data.dcd_common','sasmol.test_sasmol.data.sasmol','sasmol.test_sasmol.data.sasmol.calculate','sasmol.test_sasmol.data.sasmol.file_io','sasmol.test_sasmol.data.sasmol.file_io.test-results','sasmol.test_sasmol.data.sasmol.linear_algebra','sasmol.test_sasmol.data.sasmol.system','sasmol.test_sasmol.data.sasmol.operate','sasmol.test_sasmol.data.sasmol.properties','sasmol.test_sasmol.test_calculate','sasmol.test_sasmol.test_file_io','sasmol.test_sasmol.test_linear_algebra','sasmol.test_sasmol.test_operate','sasmol.test_sasmol.test_properties','sasmol.test_sasmol.test_subset','sasmol.extensions','sasmol.extensions.dcdio','sasmol.extensions.view','sasmol.extensions.mask','sasmol.extensions.matrix_math'],

    ext_modules=[
    Extension('sasmol._dcdio',[os.path.join('src','python','extensions','dcdio','dcdio.i'),os.path.join('src','python','extensions','dcdio','dcdio.c')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol._view_vmd',[os.path.join('src','python','extensions','view','view_vmd.i'),os.path.join('src','python','extensions','view','view_vmd.c'),os.path.join('src','python','extensions','view','imd.c'),os.path.join('src','python','extensions','view','vmdsock.c')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol._mask',[os.path.join('src','python','extensions','mask','mask.i'),os.path.join('src','python','extensions','mask','mask.c')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol.foverlap',[os.path.join('src','python','extensions','overlap','foverlap.f')],include_dirs=[NUMPY_INCLUDE]),
    Extension('sasmol.matrix_math',[os.path.join('src','python','extensions','matrix_math','matrix_math.f')],include_dirs=[NUMPY_INCLUDE])],

    data_files = [ 
        ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag.pdb')]),
        ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag_200_frames.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','1ATM.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','2AAD.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','rna-1to10.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t5.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t6.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t7.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1ATM-1to2.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1ATM.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1KP8.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','2AAD.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','rna.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1PSI.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','3AAD-2chain.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','dimcd_fixed_atoms.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1CRN.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','2AAD-1to3.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','3AAD.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','rna-1to10.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','calculate') , [os.path.join('src','python','test_sasmol','data','sasmol','calculate','1ATN.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','calculate') , [os.path.join('src','python','test_sasmol','data','sasmol','calculate','1CRN-rot.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','calculate') , [os.path.join('src','python','test_sasmol','data','sasmol','calculate','1CRN-rot-shift.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Catoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','ConflictAtoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Hatoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','MisAtoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Natoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Oatoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Otheratoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Patoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','Satoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','amino_acid_sld.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','charmm27_atoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','properties') , [os.path.join('src','python','test_sasmol','data','sasmol','properties','standard_atomic_weigh.txt')]),
   
        ( os.path.join('sasmol','test_sasmol','data','sasmol','operate') , [os.path.join('src','python','test_sasmol','data','sasmol','operate','1CRN-rot-shift.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','operate') , [os.path.join('src','python','test_sasmol','data','sasmol','operate','1CRN-rot-sub.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','operate') , [os.path.join('src','python','test_sasmol','data','sasmol','operate','1CRN-rot.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','operate') , [os.path.join('src','python','test_sasmol','data','sasmol','operate','1CRN-sub.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','linear_algebra') , [os.path.join('src','python','test_sasmol','data','sasmol','linear_algebra','1CRN-rot-shift.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','linear_algebra') , [os.path.join('src','python','test_sasmol','data','sasmol','linear_algebra','1CRN-rot.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasmol') , [os.path.join('src','python','test_sasmol','data','sasmol','sasmol','1CRN-3frames.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','1AA-NoEND.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','1ATM-1.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','1ATM-2.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-END.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-END_wrong_number_atoms.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-MODEL.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-MODEL_missing_END.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-MODEL_mix_END_noterminating.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-MODEL_wrong_number_atoms.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3-MODEL_wrongnumber_mix_END.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3_MODEL.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-1to3_MODELwrong.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-2.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','2AAD-3.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','3AAD-2chain.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','3AAD.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','nef_nohis.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','nef_nohis_1.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','file_io') , [os.path.join('src','python','test_sasmol','data','sasmol','file_io','new_package_rna.pdb')]) #,

                   ]
    )

