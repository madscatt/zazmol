from setuptools import setup, Extension
import numpy

module = Extension('dcdio_module', sources=['dcdio_module.c','dcdio.c'],
                   include_dirs=[numpy.get_include()])

#                   include_dirs=[numpy.get_include()],
#                  extra_compile_args=['-DDEBUG'])

setup(
    name='dcdio_module',
    version='1.0',
    description='Python interface for DCD file operations',
    ext_modules=[module]
)
