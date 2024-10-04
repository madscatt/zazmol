from setuptools import setup, Extension
import numpy
import os

# Ensure the 'sasmol' directory exists
if not os.path.exists('sasmol'):
    os.makedirs('sasmol')

module = Extension('sasmol._dcdio',
                   sources=['dcdio_module.c', 'dcdio.c'],
                   include_dirs=[numpy.get_include()])

setup(
    name='dcdio',
    version='1.0',
    description='Python interface for DCD file operations',
    ext_modules=[module]
)
