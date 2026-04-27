from numpy.distutils.core import setup, Extension

#matrix_math_module = Extension(
#    name='fmatrix_math',
#   sources=['matrix_math.f'],  # Use a relative path to the local folder
#   language='f77'  # Specify the Fortran language
#

matrix_math_module = Extension(
    name='fmatrix_math',
    sources=['matrix_math.f'],  # Use a relative path to the local folder
    language='f77',  # Specify the Fortran language
    libraries=['openblas'],  # Link against OpenBLAS library
    library_dirs=['/opt/homebrew/Cellar/openblas/0.3.29/lib'],  # Directory where the OpenBLAS library is located
    include_dirs=['/opt/homebrew/Cellar/openblas/0.3.29/include'],  # Directory where the OpenBLAS headers are located
)

setup(
    name='fmatrix_math',
    version='1.0',
    description='Matrix multiplication module',
    ext_modules=[matrix_math_module],
)


# python setup_fortran_mm.py build_ext --inplace