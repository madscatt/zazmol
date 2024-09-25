from setuptools import setup, Extension
import numpy

module = Extension('dcdio',
                   sources=['dcdio_module.c', 'dcdio.c'],
                   include_dirs=[numpy.get_include()])

setup(name='dcdio',
      version='1.0',
      description='DCD IO module',
      ext_modules=[module])