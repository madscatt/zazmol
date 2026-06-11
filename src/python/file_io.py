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
import os
import sys
import string
import locale
import time
import sasmol.pdb_io as pdb_io
import sasmol.mmcif_io as mmcif_io
import sasmol.dcd_io as dcd_io

#	FILE_IO
#
#	12/5/2009	--	initial coding			                    :	jc
#	12/10/2009	--	doc strings 			                    :	jc
#	01/01/2011	--	added dcdio wrappers		                :	jc
#	08/26/2016	--	split dependent classes to new files        :   jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	FILE_IO is the module that contains or inherits base classes that 
	read and write atomic information from and to the hard disk,
	and (eventually) deal with logging and general file I/O operations. 

	These classes are accessed by the SasAtm class found in
	the sasmol module.

'''

class Files(pdb_io.PDB, mmcif_io.MMCIF, dcd_io.DCD):
    '''
    Composite I/O mixin that combines PDB and DCD behavior.

    ``Files`` is inherited by ``system.Atom`` and related classes. Most public
    methods are provided by the ``PDB`` and ``DCD`` parent classes.

    Examples
    --------
    >>> import sasmol.system as system
    >>> molecule = system.Molecule()
    >>> hasattr(molecule, 'read_pdb')
    True
    '''

    def __init__(self, filename, flag):
        '''
        Placeholder initializer retained for historical API compatibility.

        Parameters
        ----------
        filename
            string : file path provided by caller

        flag
            string : mode hint provided by caller
        '''
        pass

    def _detect_structure_format(self, filename):
        '''
        Detect legacy PDB versus PDBx/mmCIF for structure input dispatch.
        '''

        lower_filename = str(filename).lower()
        if lower_filename.endswith(('.cif', '.mmcif', '.cif.gz', '.mmcif.gz')):
            return 'mmcif'
        if lower_filename.endswith(('.pdb', '.ent', '.pdb.gz')):
            return 'pdb'

        with open(filename, 'r') as infile:
            for line in infile:
                stripped = line.strip()
                if stripped == '':
                    continue
                if stripped.startswith('data_') or stripped.startswith('_atom_site.'):
                    return 'mmcif'
                record_name = line[0:6].strip()
                if record_name in ('ATOM', 'HETATM', 'HEADER', 'MODEL', 'REMARK'):
                    return 'pdb'
                break

        raise ValueError('unable to detect structure file format: ' + str(filename))

    def read_structure(self, filename, format='auto', **kwargs):
        '''
        Read a molecular structure file using format dispatch.

        ``format='auto'`` preserves legacy PDB behavior for `.pdb`/`.ent`
        inputs while routing PDBx/mmCIF files to ``read_mmcif``.
        '''

        selected_format = format
        if selected_format == 'auto':
            selected_format = self._detect_structure_format(filename)

        selected_format = selected_format.lower()
        if selected_format in ('pdb', 'legacy_pdb'):
            return self.read_pdb(filename, **kwargs)
        if selected_format in ('mmcif', 'cif', 'pdbx'):
            return self.read_mmcif(filename, **kwargs)

        raise ValueError('unsupported structure file format: ' + str(format))

    def open_file(filename):
        '''
        Placeholder helper retained for compatibility with older callers.

        Parameters
        ----------
        filename
            string : file path to open

        Returns
        -------
        None
            This method is currently a stub and performs no operation.

        Notes
        -----
        Use ``read_pdb``, ``write_pdb``, ``open_dcd_read``, and
        ``open_dcd_write`` for actual file operations.
        '''
        pass
