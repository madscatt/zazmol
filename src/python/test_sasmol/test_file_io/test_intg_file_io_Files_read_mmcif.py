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

from sasmol.test_sasmol.utilities import env

from unittest import main
import unittest
import os
import numpy
import tempfile

import sasmol.system as system


floattype = os.environ['SASMOL_FLOATTYPE']

DataRoot = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data')
PdbDataPath = os.path.join(DataRoot, 'pdb_common') + os.path.sep
MmcifDataPath = os.path.join(DataRoot, 'mmcif_common') + os.path.sep


class Test_intg_file_io_Files_read_mmcif(unittest.TestCase):
    '''
    Characterize base PDBx/mmCIF reader behavior.

    These tests intentionally cover molecule-descriptor parity only. Rich
    metadata interpretation for pdbscan/pdbrx is tracked separately.
    '''

    def assert_list_almost_equal(self, a, b, places=5):
        if len(a) != len(b):
            raise TypeError
        for i in range(len(a)):
            if isinstance(a[i], (int, float, numpy.generic)):
                if numpy.isnan(a[i]) and numpy.isnan(b[i]):
                    continue
                self.assertAlmostEqual(a[i], b[i], places)
            else:
                self.assertEqual(a[i], b[i])

    def write_temp_mmcif(self, residues):
        handle = tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False)
        try:
            handle.write('data_moltype_regression\n')
            handle.write('loop_\n')
            handle.write('_atom_site.group_PDB\n')
            handle.write('_atom_site.id\n')
            handle.write('_atom_site.type_symbol\n')
            handle.write('_atom_site.label_atom_id\n')
            handle.write('_atom_site.label_alt_id\n')
            handle.write('_atom_site.label_comp_id\n')
            handle.write('_atom_site.label_asym_id\n')
            handle.write('_atom_site.label_entity_id\n')
            handle.write('_atom_site.label_seq_id\n')
            handle.write('_atom_site.pdbx_PDB_ins_code\n')
            handle.write('_atom_site.Cartn_x\n')
            handle.write('_atom_site.Cartn_y\n')
            handle.write('_atom_site.Cartn_z\n')
            handle.write('_atom_site.occupancy\n')
            handle.write('_atom_site.B_iso_or_equiv\n')
            handle.write('_atom_site.pdbx_formal_charge\n')
            handle.write('_atom_site.auth_seq_id\n')
            handle.write('_atom_site.auth_comp_id\n')
            handle.write('_atom_site.auth_asym_id\n')
            handle.write('_atom_site.auth_atom_id\n')
            handle.write('_atom_site.pdbx_PDB_model_num\n')
            for atom_id, resname, chain in residues:
                handle.write(
                    'ATOM %d P P . %s %s 1 %d ? %.1f %.1f %.1f 1.00 0.00 ? %d %s %s P 1\n' %
                    (atom_id, resname, chain, atom_id, atom_id, atom_id + 1.0,
                     atom_id + 2.0, atom_id, resname, chain))
        finally:
            handle.close()
        return handle.name

    def test_resname_moltype_classification_matches_pdb_overlap_policy(self):
        residues = [
            (1, 'ALA', 'P'),
            (2, 'ADE', 'D'),
            (3, 'GUA', 'D'),
            (4, 'CYT', 'D'),
            (5, 'THY', 'D'),
            (6, 'URA', 'R'),
            (7, 'A', 'R'),
            (8, 'DA', 'D'),
            (9, 'TIP3', 'W'),
            (10, 'LIG', 'L')]
        cif_file = self.write_temp_mmcif(residues)
        try:
            mol = system.Molecule()
            mol.read_mmcif(cif_file)
        finally:
            os.unlink(cif_file)

        self.assertEqual(
            list(mol.moltype()),
            ['protein', 'nucleic', 'nucleic', 'nucleic', 'dna',
             'rna', 'rna', 'dna', 'water', 'other'])

    def test_1crn_mmcif_matches_existing_pdb_core_fields(self):
        pdb_mol = system.Molecule()
        cif_mol = system.Molecule()

        pdb_mol.read_pdb(PdbDataPath + '1CRN.pdb')
        cif_mol.read_mmcif(MmcifDataPath + '1CRN.cif')

        self.assertEqual(cif_mol.natoms(), pdb_mol.natoms())
        self.assertEqual(cif_mol.number_of_frames(), pdb_mol.number_of_frames())
        self.assertEqual(list(cif_mol.index()), list(pdb_mol.index()))
        self.assertEqual(list(cif_mol.original_index()), list(pdb_mol.original_index()))
        self.assertEqual(list(cif_mol.name()), list(pdb_mol.name()))
        self.assertEqual(list(cif_mol.resname()), list(pdb_mol.resname()))
        self.assertEqual(list(cif_mol.chain()), list(pdb_mol.chain()))
        self.assertEqual(list(cif_mol.resid()), list(pdb_mol.resid()))
        self.assertEqual(list(cif_mol.original_resid()), list(pdb_mol.original_resid()))
        self.assertEqual(list(cif_mol.rescode()), list(pdb_mol.rescode()))
        self.assertEqual(list(cif_mol.segname()), list(pdb_mol.segname()))
        self.assertEqual(list(cif_mol.element()), list(pdb_mol.element()))
        numpy.testing.assert_allclose(cif_mol.coor(), pdb_mol.coor(), atol=1.0e-3)

    def test_read_structure_dispatches_mmcif(self):
        mol = system.Molecule()
        mol.read_structure(MmcifDataPath + '1CRN.cif')

        self.assertEqual(mol.natoms(), 327)
        self.assertEqual(mol.number_of_frames(), 1)
        self.assertEqual(mol.name()[0], 'N')
        self.assertEqual(mol.resname()[0], 'THR')
        self.assertEqual(mol.chain()[0], 'A')

    def test_constructor_dispatches_mmcif(self):
        mol = system.Molecule(MmcifDataPath + '1CRN.cif')

        self.assertEqual(mol.natoms(), 327)
        self.assertEqual(mol.filename(), MmcifDataPath + '1CRN.cif')
        self.assertEqual(mol.number_of_frames(), 1)

    def test_nucleic_acid_and_heterogen_fixtures_read(self):
        dna = system.Molecule()
        rna = system.Molecule()
        hemoglobin = system.Molecule()

        dna.read_mmcif(MmcifDataPath + '1BNA.cif')
        rna.read_mmcif(MmcifDataPath + '1EHZ.cif')
        hemoglobin.read_mmcif(MmcifDataPath + '4HHB.cif')

        self.assertEqual(dna.natoms(), 566)
        self.assertTrue('dna' in dna.moltypes())
        self.assertEqual(rna.natoms(), 1821)
        self.assertTrue('rna' in rna.moltypes())
        self.assertEqual(hemoglobin.natoms(), 4779)
        self.assertTrue('protein' in hemoglobin.moltypes())
        self.assertTrue('HEM' in hemoglobin.resnames())

    def test_multimodel_nmr_fixture_groups_frames(self):
        mol = system.Molecule()
        mol.read_mmcif(MmcifDataPath + '2K39.cif')

        self.assertEqual(mol.natoms(), 1231)
        self.assertEqual(mol.number_of_frames(), 116)
        self.assertEqual(mol.coor().shape, (116, 1231, 3))
        self.assertEqual(list(mol.index())[:3], [1, 2, 3])


if __name__ == '__main__':
    unittest.main()
