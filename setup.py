'''
    SASMOL  Copyright (C) 2009-2016 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import numpy
import os
import subprocess
# import setuptools
# from numpy.distutils.core import Extension, setup
from setuptools import Extension, setup
from setuptools.command.build_py import build_py as _build_py
#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/05/2015      --      for sasmol distribution     :       jc
#       08/08/2016      --      branched for refactoring    :       jc
#
# LC      1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      setup.py is the script to install and/or update sasmol

	> sudo python setup.py install

'''


os.environ['NPY_DISABLE_CPU_FEATURES'] = 'AVX512FP16'

NUMPY_INCLUDE = numpy.get_include()


class build_py(_build_py):
    """Build and stage the modern native C++ SasMol library."""

    def run(self):
        _build_py.run(self)
        self._build_native_sasmol()

    def _build_native_sasmol(self):
        source_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), 'standalone_cpp_sasmol'))
        if not os.path.isdir(source_dir):
            raise RuntimeError(
                "missing modern C++ SasMol source directory: %s" % source_dir)

        build_cmd = self.get_finalized_command('build')
        build_dir = os.path.abspath(
            os.path.join(build_cmd.build_temp, 'sasmol_native'))
        stage_dir = os.path.abspath(
            os.path.join(self.build_lib, 'sasmol', 'native'))
        cmake_args = [
            'cmake',
            '-S', source_dir,
            '-B', build_dir,
            '-DCMAKE_BUILD_TYPE=Release',
            '-DBUILD_TESTING=OFF',
            '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
            '-DCMAKE_INSTALL_PREFIX=%s' % stage_dir,
        ]
        build_args = [
            'cmake',
            '--build', build_dir,
            '--target', 'install',
            '--config', 'Release',
        ]

        try:
            subprocess.check_call(cmake_args)
            subprocess.check_call(build_args)
        except OSError as exc:
            raise RuntimeError(
                "CMake is required to build the native SasMol library: %s" % exc)


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
      description=(
          "A library of methods to write software for molecular systems"),
      long_description=read('README.md'),
      classifiers=["Development Status :: 4 - Beta",
                   "License :: OSI Approved :: GNU Public License 3",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Operating System :: Linux :: MacOS :: MacOS X",
                   "Programming Language :: Python :: C :: Fortran",
                   "Topic :: Scientific/Engineering :: Chemistry :: Physics"],

      package_dir={
          'sasmol': os.path.join('src', 'python'),
          'mmcif': os.path.join(
              'third_party', 'rcsb', 'mmcif-1.1.1', 'mmcif')},

      packages=['sasmol', 'mmcif', 'mmcif.api', 'mmcif.core', 'mmcif.io', 'sasmol.test_sasmol', 'sasmol.test_sasmol.utilities', 'sasmol.test_sasmol.data', 'sasmol.test_sasmol.data.pdb_common', 'sasmol.test_sasmol.data.dcd_common', 'sasmol.test_sasmol.data.sasmol', 'sasmol.test_sasmol.data.sasmol.calculate', 'sasmol.test_sasmol.data.sasmol.file_io', 'sasmol.test_sasmol.data.sasmol.file_io.test-results', 'sasmol.test_sasmol.data.sasmol.linear_algebra', 'sasmol.test_sasmol.test_utilities', 'sasmol.test_sasmol.test_system', 'sasmol.test_sasmol.data.sasmol.system',
                'sasmol.test_sasmol.data.sasmol.operate', 'sasmol.test_sasmol.data.sasmol.properties', 'sasmol.test_sasmol.test_calculate', 'sasmol.test_sasmol.test_file_io', 'sasmol.test_sasmol.test_linear_algebra', 'sasmol.test_sasmol.test_operate', 'sasmol.test_sasmol.test_properties', 'sasmol.test_sasmol.test_subset', 'sasmol.test_sasmol.test_topology', 'sasmol.extensions', 'sasmol.extensions.dcdio', 'sasmol.extensions.view', 'sasmol.extensions.mask', 'sasmol.extensions.matrix_math'],

      package_data={
          'sasmol': [
              'native/include/sasmol/*.hpp',
              'native/lib/*.a',
              'native/lib/cmake/sasmol/*.cmake',
          ],
      },

      cmdclass={'build_py': build_py},

      ext_modules=[
          Extension('sasmol._dcdio',
                    sources=[os.path.join('src', 'python', 'extensions', 'dcdio', 'dcdio.c'),
                             os.path.join('src', 'python', 'extensions', 'dcdio', 'dcdio_module.c')],
                    include_dirs=[numpy.get_include()]),

          Extension('sasmol._view_vmd', [os.path.join('src', 'python', 'extensions', 'view', 'view_vmd_module.c'), os.path.join('src', 'python', 'extensions', 'view', 'view_vmd.c'), os.path.join(
              'src', 'python', 'extensions', 'view', 'imd.c'), os.path.join('src', 'python', 'extensions', 'view', 'vmdsock.c')], include_dirs=[NUMPY_INCLUDE]),
          Extension('sasmol._mask', [os.path.join('src', 'python', 'extensions', 'mask', 'mask_module.c'), os.path.join(
              'src', 'python', 'extensions', 'mask', 'mask.c')], include_dirs=[NUMPY_INCLUDE]),
          Extension('sasmol.overlap', [os.path.join(
              'src', 'python', 'extensions', 'overlap', 'overlap.c')], include_dirs=[NUMPY_INCLUDE]),
          Extension('sasmol.matrix_math', [os.path.join('src', 'python', 'extensions', 'matrix_math', 'matrix_math.c')], include_dirs=[NUMPY_INCLUDE])],

      data_files=[
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'dcd_common', '1ATM.dcd')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'dcd_common', '2AAD.dcd')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'dcd_common', 'rna-1to10.dcd')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [
           os.path.join('src', 'python', 'test_sasmol', 'data', 'dcd_common', 't.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [
           os.path.join('src', 'python', 'test_sasmol', 'data', 'dcd_common', 't5.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [
           os.path.join('src', 'python', 'test_sasmol', 'data', 'dcd_common', 't6.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'dcd_common'), [
           os.path.join('src', 'python', 'test_sasmol', 'data', 'dcd_common', 't7.pdb')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '1ATM-1to2.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '1ATM.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '1KP8.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '2AAD.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [
              os.path.join('src', 'python', 'test_sasmol', 'data', 'pdb_common', 'rna.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '1PSI.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '3AAD-2chain.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', 'dimcd_fixed_atoms.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '1CRN.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '2AAD-1to3.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', '3AAD.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'pdb_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'pdb_common', 'rna-1to10.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', '1BNA.cif')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', '1CRN.cif')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', '1EHZ.cif')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', '1KP8.cif')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', '2K39.cif')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', '4HHB.cif')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'mmcif_common'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'mmcif_common', 'README.md')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'calculate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'calculate', '1ATN.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'calculate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'calculate', '1CRN-rot.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'calculate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'calculate', '1CRN-rot-shift.pdb')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Catoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'ConflictAtoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Hatoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'MisAtoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Natoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Oatoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Otheratoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Patoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'Satoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'amino_acid_sld.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'charmm27_atoms.txt')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'properties'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'properties', 'standard_atomic_weight.txt')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'operate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'operate', '1CRN-rot-shift.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'operate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'operate', '1CRN-rot-sub.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'operate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'operate', '1CRN-rot.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'operate'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'operate', '1CRN-sub.pdb')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'linear_algebra'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'linear_algebra', '1CRN-rot-shift.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'linear_algebra'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'linear_algebra', '1CRN-rot.pdb')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'sasmol'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'system', '1CRN-3frames.pdb')]),

          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '1AA-NoEND.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '1ATM-1.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '1ATM-2.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-END.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-END_wrong_number_atoms.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-MODEL.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-MODEL_missing_END.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-MODEL_mix_END_noterminating.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-MODEL_wrong_number_atoms.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3-MODEL_wrongnumber_mix_END.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3_MODEL.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-1to3_MODELwrong.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-2.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '2AAD-3.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '3AAD-2chain.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', '3AAD.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', 'nef_nohis.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', 'nef_nohis_1.pdb')]),
          (os.path.join('sasmol', 'test_sasmol', 'data', 'sasmol', 'file_io'), [os.path.join(
              'src', 'python', 'test_sasmol', 'data', 'sasmol', 'file_io', 'new_package_rna.pdb')])  # ,

      ]
      )
