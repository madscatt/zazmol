'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# simple extension module
_mask = Extension("_mask",
                   ["mask.i","mask.c"],
                   include_dirs = [numpy_include],
                   swig_opts=['-py3'],
                   )

# NumyTypemapTests setup
setup(  name        = "MASK",
        description = "Module handles mask selection",
        author      = "Joseph E. Curtis",
        version     = "0.1",
        ext_modules = [_mask]
        )

