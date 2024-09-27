from setuptools import setup, Extension

module = Extension('dcdio_module', sources=['dcdio_module.c','dcdio.c'])

setup(
    name='dcdio_module',
    version='1.0',
    description='Python interface for DCD file operations',
    ext_modules=[module]
)
