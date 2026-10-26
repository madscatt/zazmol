"""Reader metadata contracts; no chemical interpretation is performed."""
import os
import tempfile
import unittest
from unittest import mock

import numpy
import sasmol.config as config
import sasmol.system as system


ATTRIBUTES = '''group_PDB id type_symbol label_atom_id auth_atom_id label_alt_id
label_comp_id auth_comp_id label_asym_id auth_asym_id label_entity_id
label_seq_id auth_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z
occupancy B_iso_or_equiv pdbx_PDB_model_num'''.split()
ROWS = [
    'ATOM 1 N N N . ALA ALA L1 A 1 1 10 A 0 1 2 1 20 7'.split(),
    'ATOM 1 C CA CA ? ALA ALA L2 A 2 1 10 A 3 4 5 ? 20 7'.split(),
    'ATOM 1 C CA CA B ALA ALA L2 A 2 1 10 A 6 7 8 . 20 9'.split(),
    'ATOM 1 N N N . ALA ALA L1 A 1 1 10 A 9 10 11 0 20 9'.split(),
]
EXTRA = '''
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.pdbx_ptnr2_label_alt_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
c1 covale L1 L2 B 1_555 2_565
#
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
#
_struct_asym.id L1
_struct_asym.entity_id 1
_entity_poly.entity_id 1
_entity_poly.type 'polypeptide(L)'
_entity_poly_seq.entity_id 1
_entity_poly_seq.num 1
_entity_poly_seq.mon_id ALA
_pdbx_poly_seq_scheme.asym_id L1
_pdbx_poly_seq_scheme.seq_id 1
_pdbx_poly_seq_scheme.auth_seq_num 10
_pdbx_nonpoly_scheme.asym_id L2
_pdbx_nonpoly_scheme.auth_seq_num 10
_pdbx_unobs_or_zero_occ_residues.polymer_flag Y
_pdbx_unobs_or_zero_occ_residues.label_seq_id 2
_pdbx_unobs_or_zero_occ_atoms.polymer_flag N
_pdbx_unobs_or_zero_occ_atoms.label_seq_id .
_pdbx_unobs_or_zero_occ_atoms.occupancy_flag 0
_pdbx_struct_assembly.id 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1-2)(3,4)' L1
1 1 L2
#
_pdbx_struct_oper_list.id 1
_pdbx_struct_oper_list.matrix[1][1] 1.0
_custom_evidence.value 'uninterpreted text'
data_other
_entry.id SECOND
'''


def cif_text(rows=None, omit=()):
    attributes = [a for a in ATTRIBUTES if a not in omit]
    indices = [ATTRIBUTES.index(a) for a in attributes]
    return ('data_probe\nloop_\n' +
            ''.join('_atom_site.' + a + '\n' for a in attributes) +
            '\n'.join(' '.join(row[i] for i in indices)
                      for row in (ROWS if rows is None else rows)) + '\n#\n')


class Test_mmcif_metadata(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = self.write('source.cif', cif_text() + EXTRA)

    def write(self, name, text):
        path = os.path.join(self.directory.name, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def read(self, path=None, **kwargs):
        mol = system.Molecule()
        self.assertIsNone(mol.read_mmcif(path or self.path, **kwargs))
        return mol

    def test_raw_categories_namespaces_nulls_and_frame_rows(self):
        mol = self.read(pdbscan=True)
        metadata = mol.mmcif_metadata()
        self.assertEqual(metadata['schema_version'], 1)
        self.assertEqual(metadata['filename'], self.path)
        self.assertEqual(metadata['coordinate_block_index'], 0)
        self.assertEqual(metadata['correspondence_status'], 'read_order')
        self.assertEqual(metadata['frame_atom_site_rows'], [[0, 1], [2, 3]])
        self.assertEqual([b['name'] for b in metadata['data_blocks']],
                         ['probe', 'other'])
        cats = metadata['data_blocks'][0]['categories']
        self.assertEqual(cats['atom_site']['attributes'], ATTRIBUTES)
        self.assertEqual(cats['atom_site']['rows'], ROWS)
        self.assertEqual(cats['struct_conn']['rows'],
                         [['c1', 'covale', 'L1', 'L2', 'B', '1_555', '2_565']])
        self.assertEqual(cats['pdbx_struct_assembly_gen']['rows'],
                         [['1', '(1-2)(3,4)', 'L1'], ['1', '1', 'L2']])
        self.assertEqual(cats['pdbx_unobs_or_zero_occ_atoms']['rows'],
                         [['N', '.', '0']])
        for category in ('entity', 'struct_asym', 'entity_poly', 'entity_poly_seq',
                         'pdbx_poly_seq_scheme', 'pdbx_nonpoly_scheme',
                         'pdbx_unobs_or_zero_occ_residues', 'pdbx_struct_assembly',
                         'pdbx_struct_oper_list', 'custom_evidence'):
            self.assertIn(category, cats)
        self.assertEqual(metadata['data_blocks'][1]['categories']['entry']['rows'],
                         [['SECOND']])

    def test_absent_fields_are_not_fabricated(self):
        path = self.write('absent.cif', cif_text(omit=(
            'label_alt_id', 'occupancy', 'pdbx_PDB_model_num')))
        metadata = self.read(path).mmcif_metadata()
        attrs = metadata['data_blocks'][0]['categories']['atom_site']['attributes']
        for attr in ('label_alt_id', 'occupancy', 'pdbx_PDB_model_num'):
            self.assertNotIn(attr, attrs)
        self.assertEqual(metadata['frame_atom_site_rows'], [[0, 1, 2, 3]])

    def test_existing_coordinates_descriptors_and_scan_defaults(self):
        for scan in (False, True):
            with self.subTest(pdbscan=scan):
                mol = self.read(pdbscan=scan)
                numpy.testing.assert_array_equal(mol.coor(),
                    numpy.arange(12).reshape(2, 2, 3))
                self.assertEqual(mol.coor().dtype, numpy.dtype(config.COORD_DTYPE))
                self.assertEqual(mol.name(), ['N', 'CA'])
                self.assertEqual(mol.chain(), ['A', 'A'])
                self.assertEqual(mol.segname(), ['A', 'A'])
                numpy.testing.assert_array_equal(mol.resid(), [10, 10])
                numpy.testing.assert_array_equal(mol.original_index(), [1, 1])
                numpy.testing.assert_array_equal(mol.index(), [1, 2])
                self.assertEqual(mol.loc(), [' ', ' '])
                self.assertEqual(mol.occupancy(), ['1', '' if scan else '  1.00'])
                self.assertEqual(mol.header(), [])
                self.assertEqual(mol.conect(), {})
                self.assertEqual(mol.biomt(), {})
        rows = [row[:] for row in ROWS[:2]]
        rows[1][5] = 'B'
        path = self.write('alternate.cif', cif_text(rows))
        self.assertEqual(self.read(path, pdbscan=True).loc(), [' ', 'B'])
        self.assertEqual(self.read(path).loc(), [' ', ' '])

    def test_snapshot_isolated_from_accessor_and_working_edits(self):
        mol = self.read()
        expected = mol.mmcif_metadata()
        result = mol.mmcif_metadata()
        result['data_blocks'][0]['categories']['atom_site']['rows'][0][0] = 'EDIT'
        result['frame_atom_site_rows'][0].reverse()
        mol.setChain(['X', 'Y'])
        mol.setSegname(['XX', 'YY'])
        mol.coor()[:] += 100
        self.assertEqual(mol.mmcif_metadata(), expected)
        duplicate = mol.duplicate_molecule(1)[0]
        self.assertEqual(duplicate.mmcif_metadata(), expected)
        mol.read_mmcif(self.write('replacement.cif', cif_text(ROWS[:1])))
        self.assertEqual(duplicate.mmcif_metadata(), expected)
        self.assertEqual(mol.mmcif_metadata()['frame_atom_site_rows'], [[0]])

    def test_fresh_pdb_and_failed_reads_have_no_snapshot(self):
        self.assertIsNone(system.Molecule().mmcif_metadata())
        mol = self.read()
        # Use the established PDB fixture for successful reuse.
        pdb_path = os.path.join(os.path.dirname(__file__), '..', 'data',
                                'pdb_common', '1CRN.pdb')
        mol.read_pdb(pdb_path)
        self.assertIsNone(mol.mmcif_metadata())
        for suffix, reader in (('cif', 'read_mmcif'), ('pdb', 'read_pdb')):
            mol = self.read()
            with self.assertRaises(Exception):
                getattr(mol, reader)(os.path.join(self.directory.name, 'missing.' + suffix))
            self.assertIsNone(mol.mmcif_metadata())
        mol = self.read()
        with self.assertRaises(ValueError):
            mol.read_mmcif(self.write('unequal.cif', cif_text(ROWS[:3])))
        self.assertIsNone(mol.mmcif_metadata())

    def test_first_block_and_topology_behavior(self):
        mol = self.read()
        with self.assertRaises(ValueError):
            mol.read_mmcif(self.write('blocks.cif', 'data_empty\n_entry.id EMPTY\n' + cif_text()))
        self.assertIsNone(mol.mmcif_metadata())
        for result in ([], ['existing topology diagnostic']):
            with mock.patch.object(mol, 'check_charmm_atomic_order_reorganize',
                                   return_value=result) as topology:
                self.assertIs(mol.read_mmcif(self.path, saspdbrx_topology=True), result)
                topology.assert_called_once_with()
            metadata = mol.mmcif_metadata()
            self.assertEqual(metadata['correspondence_status'], 'unavailable_after_topology')
            self.assertIsNone(metadata['frame_atom_site_rows'])
            self.assertEqual(metadata['data_blocks'][0]['categories']['atom_site']['rows'], ROWS)
        with mock.patch.object(mol, 'check_charmm_atomic_order_reorganize',
                               side_effect=ValueError('existing failure')):
            with self.assertRaises(ValueError):
                mol.read_mmcif(self.path, saspdbrx_topology=True)
        self.assertIsNone(mol.mmcif_metadata())

    def test_subset_and_merge_clear_only_destination_metadata(self):
        source, target = self.read(), self.read()
        expected = source.mmcif_metadata()
        self.assertEqual(source.copy_molecule_using_mask(target, [1, 0], 1), [])
        self.assertIsNone(target.mmcif_metadata())
        numpy.testing.assert_array_equal(target.coor(), [[[6, 7, 8]]])
        self.assertEqual(source.mmcif_metadata(), expected)
        target = self.read()
        self.assertEqual(target.merge_two_molecules(source, source), [])
        self.assertIsNone(target.mmcif_metadata())
        self.assertEqual(source.mmcif_metadata(), expected)
        numpy.testing.assert_array_equal(target.coor(), [[[0, 1, 2], [3, 4, 5],
                                                         [0, 1, 2], [3, 4, 5]]])


if __name__ == '__main__':
    unittest.main()
