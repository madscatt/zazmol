# from __future__ import unicode_literals
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

import gzip
import os
import sys

import numpy

from . import config as config
from .pdb_io import PDBElementResolutionError


class MMCIFIdentifierError(ValueError):
    '''Raised when mmCIF identifiers cannot fit SASMOL descriptor semantics.'''


def _ensure_mmcif_parser_available():
    '''
    Make the official RCSB/wwPDB parser snapshot importable during source-tree
    development while still allowing an installed ``mmcif`` dependency.
    '''

    try:
        from mmcif.io.PdbxReader import PdbxReader
        return PdbxReader
    except ImportError:
        pass

    source_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    vendored_path = os.path.join(source_root, 'third_party', 'rcsb', 'mmcif-1.1.1')
    if vendored_path not in sys.path and os.path.isdir(vendored_path):
        sys.path.insert(0, vendored_path)

    from mmcif.io.PdbxReader import PdbxReader
    return PdbxReader


def _is_missing_mmcif_value(value):
    return value is None or value == '' or value == '.' or value == '?'


def _clean_mmcif_value(value, default=''):
    if _is_missing_mmcif_value(value):
        return default
    return value


class MMCIF(object):
    '''
    Methods for reading PDBx/mmCIF coordinate data into SASMOL descriptors.

    The reader intentionally targets the same base molecule state as
    ``read_pdb()``. Rich metadata interpretation for ``pdbscan``/``pdbrx`` is a
    separate workflow and should not be folded into this base reader.
    '''

    def _read_mmcif_data_blocks(self, filename):
        PdbxReader = _ensure_mmcif_parser_available()
        data = []

        if filename.endswith('.gz'):
            with gzip.open(filename, 'rt') as infile:
                PdbxReader(infile).read(data)
        else:
            with open(filename, 'r') as infile:
                PdbxReader(infile).read(data)

        if len(data) == 0:
            raise ValueError('mmCIF file did not contain a data block: ' + filename)

        return data

    def _mmcif_value(self, category, attribute_name, row_index, default=''):
        if category is None or not category.hasAttribute(attribute_name):
            return default
        return _clean_mmcif_value(category.getValue(attribute_name, row_index), default)

    def _mmcif_first_value(self, category, attribute_names, row_index, default=''):
        for attribute_name in attribute_names:
            value = self._mmcif_value(category, attribute_name, row_index, None)
            if not _is_missing_mmcif_value(value):
                return value
        return default

    def _mmcif_int_value(self, category, attribute_names, row_index, label):
        value = self._mmcif_first_value(category, attribute_names, row_index, None)
        if _is_missing_mmcif_value(value):
            raise MMCIFIdentifierError('missing mmCIF ' + label)

        try:
            return int(value)
        except Exception:
            raise MMCIFIdentifierError(
                'mmCIF ' + label + ' is not compatible with SASMOL integer descriptors: ' +
                str(value))

    def _mmcif_model_value(self, atom_site, row_index):
        return self._mmcif_value(atom_site, 'pdbx_PDB_model_num', row_index, '1')

    def _mmcif_model_groups(self, atom_site):
        model_order = []
        model_rows = {}

        for row_index in range(atom_site.getRowCount()):
            model = self._mmcif_model_value(atom_site, row_index)
            if model not in model_rows:
                model_order.append(model)
                model_rows[model] = []
            model_rows[model].append(row_index)

        return [model_rows[model] for model in model_order]

    def _mmcif_required_atom_site(self, data_block):
        atom_site = data_block.getObj('atom_site')
        if atom_site is None or atom_site.getRowCount() == 0:
            raise ValueError('mmCIF file does not contain _atom_site records')
        return atom_site

    def _mmcif_moltype_for_resname(self, resname):
        protein_resnames, dna_resnames, rna_resnames, nucleic_resnames, water_resnames = self.get_resnames()

        if resname in protein_resnames:
            return 'protein'
        if resname in rna_resnames:
            return 'rna'
        if resname in dna_resnames:
            return 'dna'
        if resname in water_resnames:
            return 'water'
        return 'other'

    def _mmcif_unique_counts(self, values):
        unique_values = []
        for value in values:
            if value not in unique_values:
                unique_values.append(value)
        return unique_values

    def read_mmcif(self, filename, **kwargs):
        '''
        Read a PDBx/mmCIF file and populate molecule descriptors.

        This base reader maps ``_atom_site`` coordinate data into the same
        operational descriptors as ``read_pdb()``. It deliberately does not
        interpret ``_struct_conn`` or other metadata categories for
        ``pdbscan``/``pdbrx``.
        '''

        data = self._read_mmcif_data_blocks(filename)
        atom_site = self._mmcif_required_atom_site(data[0])
        model_groups = self._mmcif_model_groups(atom_site)

        num_frames = len(model_groups)
        num_atoms = len(model_groups[0])
        for rows in model_groups:
            if len(rows) != num_atoms:
                raise ValueError('number of atoms per mmCIF model is not equal')

        coor = numpy.zeros((num_frames, num_atoms, 3), config.COORD_DTYPE)

        first_model = model_groups[0]
        atom = []
        index = []
        original_index = []
        name = []
        loc = []
        resname = []
        chain = []
        resid = []
        original_resid = []
        rescode = []
        occupancy = []
        beta = []
        segname = []
        element = []
        charge = []
        moltype = []
        residue_flag = []

        pdbscan = kwargs.get('pdbscan', False)

        for atom_ordinal, row_index in enumerate(first_model):
            record_name = self._mmcif_value(atom_site, 'group_PDB', row_index, 'ATOM').strip()
            atom.append(record_name)

            original_index.append(
                self._mmcif_int_value(atom_site, ['id'], row_index, 'atom id'))
            index.append(atom_ordinal + 1)

            this_name = self._mmcif_first_value(
                atom_site, ['auth_atom_id', 'label_atom_id'], row_index, '').strip()
            name.append(this_name)

            altloc = self._mmcif_value(atom_site, 'label_alt_id', row_index, ' ')
            loc.append(altloc if pdbscan else ' ')

            this_resname = self._mmcif_first_value(
                atom_site, ['auth_comp_id', 'label_comp_id'], row_index, '').strip()
            resname.append(this_resname)

            this_chain = self._mmcif_first_value(
                atom_site, ['auth_asym_id', 'label_asym_id'], row_index, '').strip()
            chain.append(this_chain)

            this_resid = self._mmcif_int_value(
                atom_site, ['auth_seq_id', 'label_seq_id'], row_index, 'residue id')
            resid.append(this_resid)
            original_resid.append(this_resid)

            this_rescode = self._mmcif_value(atom_site, 'pdbx_PDB_ins_code', row_index, ' ')
            rescode.append(this_rescode)

            this_occupancy = self._mmcif_value(atom_site, 'occupancy', row_index, None)
            if _is_missing_mmcif_value(this_occupancy):
                this_occupancy = '' if pdbscan else '  1.00'
            occupancy.append(this_occupancy)

            this_beta = self._mmcif_value(atom_site, 'B_iso_or_equiv', row_index, None)
            if _is_missing_mmcif_value(this_beta):
                this_beta = '' if pdbscan else '  0.00'
            beta.append(this_beta)

            # mmCIF has no direct legacy PDB segment-id counterpart.
            segname.append(this_chain)

            this_element = self._mmcif_value(atom_site, 'type_symbol', row_index, '  ').strip()
            if not this_element:
                this_element = '  '
            element.append(this_element.upper())

            this_charge = self._mmcif_value(atom_site, 'pdbx_formal_charge', row_index, '  ')
            if _is_missing_mmcif_value(this_charge):
                this_charge = '  '
            charge.append(this_charge)

            moltype.append(self._mmcif_moltype_for_resname(this_resname))
            residue_flag.append(False)

        for frame_index, rows in enumerate(model_groups):
            for atom_index, row_index in enumerate(rows):
                try:
                    coor[frame_index, atom_index, 0] = config.COORD_DTYPE(
                        self._mmcif_value(atom_site, 'Cartn_x', row_index))
                    coor[frame_index, atom_index, 1] = config.COORD_DTYPE(
                        self._mmcif_value(atom_site, 'Cartn_y', row_index))
                    coor[frame_index, atom_index, 2] = config.COORD_DTYPE(
                        self._mmcif_value(atom_site, 'Cartn_z', row_index))
                except Exception:
                    raise ValueError('failed to parse mmCIF coordinates')

        self._atom = atom
        self._index = numpy.array(index, int)
        self._original_index = numpy.array(original_index, int)
        self._name = name
        self._loc = loc
        self._resname = resname
        self._chain = chain
        self._resid = numpy.array(resid, int)
        self._original_resid = numpy.array(original_resid, int)
        self._rescode = rescode
        self._occupancy = occupancy
        self._beta = beta
        self._segname = segname
        self._element = element
        self._charge = charge
        self._moltype = moltype
        self._residue_flag = residue_flag
        self._coor = numpy.array(coor)

        if 'check_zero_coor' in kwargs:
            self.check_for_all_zero_columns(self._coor)

        try:
            unique_elements = self.element_filter()
        except PDBElementResolutionError:
            raise

        unique_names = self._mmcif_unique_counts(name)
        unique_resnames = self._mmcif_unique_counts(resname)
        unique_resids = self._mmcif_unique_counts(resid)
        unique_chains = self._mmcif_unique_counts(chain)
        unique_segnames = self._mmcif_unique_counts(segname)
        unique_occupancies = self._mmcif_unique_counts(occupancy)
        unique_betas = self._mmcif_unique_counts(beta)
        unique_moltypes = self._mmcif_unique_counts(moltype)

        self._number_of_names = len(unique_names)
        self._names = unique_names
        self._number_of_resnames = len(unique_resnames)
        self._resnames = unique_resnames
        self._number_of_resids = len(unique_resids)
        self._resids = unique_resids
        self._number_of_chains = len(unique_chains)
        self._chains = unique_chains
        self._number_of_segnames = len(unique_segnames)
        self._segnames = unique_segnames
        self._number_of_occupancies = len(unique_occupancies)
        self._occupancies = unique_occupancies
        self._number_of_betas = len(unique_betas)
        self._betas = unique_betas
        self._number_of_moltypes = len(unique_moltypes)
        self._moltypes = unique_moltypes
        self._number_of_elements = len(unique_elements)
        self._elements = unique_elements
        self._natoms = len(index)

        self.calculate_mass()

        if 'saspdbrx_topology' in kwargs:
            if kwargs['saspdbrx_topology']:
                return self.check_charmm_atomic_order_reorganize()

        self._header = []
        self._biomt = {}
        self._conect = {}

        return
