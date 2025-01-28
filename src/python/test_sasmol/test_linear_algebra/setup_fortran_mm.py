from setuptools import setup, Extension
from numpy.distutils.core import setup, Extension

matrix_math_module = Extension(
    name='matrix_math',
    sources=['/Users/curtisj/git_working_copies/zazmol/src/python/test_sasmol/test_linear_algebra/matrix_math.f'],
)

setup(
    name='matrix_math',
    version='1.0',
    description='Matrix multiplication module',
    ext_modules=[matrix_math_module],
)