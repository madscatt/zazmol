from setuptools import setup, Extension

overlap = Extension('overlap',
                    sources = ['overlap.c'])

setup(  name        = "OVERLAP",
        description = "Module determines overlap between atoms",
        author      = "Joseph E. Curtis",
        version     = "1.0",
        ext_modules = [overlap]
        )
