from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
#
'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
#### NOTE: INDENTATION FIXED : 12/11/2015 JEC

import os
import sys
import string
import locale
import struct
import numpy
import time
import sasmol.dcdio as dcdio
import sasmol.pdb_io as pdb_io
import sasmol.dcd_io as dcd_io

#	FILE_IO
#
#	12/5/2009	--	initial coding			:	jc
#	12/10/2009	--	doc strings 			:	jc
#	01/01/2011	--	added dcdio wrappers		:	jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	FILE_IO is the main module that contains the base classes that 
	read and write atomic information from and to the hard disk,
	and (eventually) deal with logging of input parameters and
	project runs.

	The methods in class Files are used to read and write data to
	the Charmm/Xplor binary data format (DCD) and textual protein
	data bank (PDB) format.  

	See the following sites for the DCD format:

	http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
	http://www.lrz-muenchen.de/~heller/ego/manual/node93.html

	The methods (read_pdb, write_pdb) are used to read and write data
	using the PDB file format as described at the following
	web-site:

	http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

	These classes are accessed by the SasAtm class found in
	the sasmol module.

'''

class Files(pdb_io.PDB, dcd_io.DCD):

    def __init__(self,filename,flag):
        pass


