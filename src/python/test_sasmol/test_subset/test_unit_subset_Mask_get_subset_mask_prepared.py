"""Differential checks for the separate prepared-selection API (issue #106)."""
import inspect
import pickle
from pathlib import Path
import unittest

import numpy
from sasmol import subset, system


class Test_prepared_subset_mask(unittest.TestCase):
    def setUp(self):
        self.molecule = system.Molecule_Maker(
            6, index=[1, 2, 3, 4, 5, 6],
            name=['N', 'CA', 'P', 'O', 'C1', 'P'],
            loc=[' ', 'A', ' ', 'B', ' ', ' '],
            resname=['ALA', 'ALA', 'DA', 'HOH', 'LIG', 'U'],
            chain=['A', 'A', 'B', 'A', 'C', 'C'],
            resid=[9, 9, -1, 9, 2, 2],
            rescode=['', 'A', '', '', 'B', ''],
            occupancy=['1.00', '0.50', '1.00', '1.00', '0.00', '1.00'],
            beta=['10.0', '20.0', '0.0', '0.0', '0.0', '0.0'],
            segname=['PRO', 'PRO', 'DNA', '', 'LIG', 'RNA'],
            element=['N', 'C', 'P', 'O', 'C', 'P'],
            charge=['', '', '', '', '1+', ''],
            moltype=['protein', 'protein', 'dna', 'water', 'other', 'rna'],
            residue_flag=[0, 1, 0, 1, 0, 1])

    def compare(self, expression, molecule=None):
        molecule = self.molecule if molecule is None else molecule
        before = pickle.dumps(vars(molecule), protocol=4)
        old_error, old_mask = molecule.get_subset_mask(expression)
        new_error, new_mask = molecule.get_subset_mask_prepared(expression)
        self.assertEqual(new_error, old_error)
        self.assertIs(type(new_mask), type(old_mask))
        numpy.testing.assert_array_equal(new_mask, old_mask)
        if isinstance(old_mask, numpy.ndarray):
            self.assertEqual(new_mask.dtype, old_mask.dtype)
            numpy.testing.assert_array_equal(
                molecule.get_indices_from_mask(new_mask),
                molecule.get_indices_from_mask(old_mask))
        self.assertEqual(pickle.dumps(vars(molecule), protocol=4), before)
        return new_error, new_mask

    def test_all_descriptors_and_combined_selections(self):
        expressions = [
            'index[i] == 2', 'name[i] == "CA"', 'loc[i] == "A"',
            'resname[i] == "ALA"', 'chain[i] == "C"', 'resid[i] == 9',
            'rescode[i] == "B"', 'occupancy[i] == "0.50"',
            'beta[i] == "20.0"', 'segname[i] == ""',
            'element[i] == "P"', 'charge[i] == "1+"',
            'moltype[i] == "other"', 'residue_flag[i] == 1',
            'chain[i] == "A" and segname[i] == "PRO" and resid[i] == 9',
            '(resid[i] == 9 and rescode[i] == "A") or moltype[i] == "rna"',
        ]
        for expression in expressions:
            with self.subTest(expression=expression):
                error, mask = self.compare(expression)
                self.assertEqual(error, [])

    def test_live_descriptors_and_successive_queries(self):
        expression = 'moltype[i] == "other"'
        self.compare(expression)
        self.molecule.setMoltype(['other'] * 6)
        error, mask = self.compare(expression)
        self.assertEqual(error, [])
        numpy.testing.assert_array_equal(mask, [1] * 6)
        self.molecule.residue_flag()[:] = [1] * 6
        self.compare('residue_flag[i] == 1')
        self.molecule.setSegname(['NEW'] * 6)
        self.compare('segname[i] == "NEW"')
        self.compare('name[i] == "P"')
        self.compare(expression)

    def test_changed_atom_count_is_read_on_next_call(self):
        self.compare('True')
        # Shrink an otherwise unchanged molecule through existing public setters.
        self.molecule.setNatoms(3)
        self.molecule.setCoor(self.molecule.coor()[:, :3, :].copy())
        for getter, setter in [
                ('index', 'setIndex'), ('name', 'setName'), ('loc', 'setLoc'),
                ('resname', 'setResname'), ('chain', 'setChain'),
                ('resid', 'setResid'), ('rescode', 'setRescode'),
                ('occupancy', 'setOccupancy'), ('beta', 'setBeta'),
                ('segname', 'setSegname'), ('element', 'setElement'),
                ('charge', 'setCharge'), ('moltype', 'setMoltype')]:
            getattr(self.molecule, setter)(getattr(self.molecule, getter)()[:3])
        self.molecule.residue_flag()[:] = self.molecule.residue_flag()[:3]
        error, mask = self.compare('True')
        self.assertEqual(error, [])
        self.assertEqual(len(mask), 3)

    def test_existing_empty_no_match_and_failure_boundaries(self):
        for expression in ['False', 'invalid =', 'i < 1 or 1 / (i - 1)']:
            with self.subTest(expression=expression):
                self.compare(expression)
        self.compare('invalid =', system.Molecule_Maker(0))

    def test_missing_descriptor_keeps_existing_exception_boundary(self):
        del self.molecule._resname
        for method in [self.molecule.get_subset_mask,
                       self.molecule.get_subset_mask_prepared]:
            with self.assertRaises(AttributeError):
                method('True')

    def test_whitespace_and_label_expression_text(self):
        for label in ['A\\n', 'A\\t', 'A"', "A'", 'A\\', 'Å']:
            self.molecule.setChain([label] * 6)
            expression = 'chain[i] == "' + label + '"'
            with self.subTest(label=label):
                self.compare(expression)
        for expression in [' \tname[i] == "N"', '\nname[i] == "N"\n']:
            self.compare(expression)

    def test_string_subclasses_use_stored_expression_text(self):
        class Expression(str):
            def __str__(self):
                return 'False'

            def lstrip(self, *args):
                return 'False'

        self.compare(Expression(' \tname[i] == "CA"'))

    def test_original_signature_and_expression_namespace_are_unchanged(self):
        self.assertEqual(str(inspect.signature(subset.Mask.get_subset_mask)),
                         '(self, basis_filter)')
        error, mask = self.molecule.get_subset_mask('prepared')
        self.assertEqual(error, ['failed to evaluate filter selection prepared for atom 0'])
        self.assertEqual(mask, [])

    def test_existing_subclass_override_is_not_redirected(self):
        class CustomMolecule(system.Molecule_Maker):
            def get_subset_mask(self, basis_filter):
                return ['custom selector'], [7]

        molecule = CustomMolecule(2, name=['N', 'CA'])
        self.assertEqual(molecule.get_subset_mask('anything'),
                         (['custom selector'], [7]))
        error, mask = molecule.get_subset_mask_prepared('name[i] == "CA"')
        self.assertEqual(error, [])
        numpy.testing.assert_array_equal(mask, [0, 1])

    def test_molecule_inherits_the_new_method(self):
        self.assertIs(system.Molecule.get_subset_mask_prepared,
                      subset.Mask.get_subset_mask_prepared)
        self.assertIs(system.Molecule_Maker.get_subset_mask_prepared,
                      subset.Mask.get_subset_mask_prepared)

    def test_reader_populated_descriptors(self):
        fixtures = Path(__file__).resolve().parents[1] / 'data' / 'pdb_common'
        for filename in ['1CRN.pdb', '1KP8.pdb', 'rna.pdb']:
            molecule = system.Molecule(0)
            molecule.read_pdb(str(fixtures / filename))
            for expression in ['name[i] == "CA"', 'chain[i] == "A"',
                               'resid[i] < 10', 'rescode[i] != ""',
                               'moltype[i] == "protein"', 'segname[i] == ""']:
                with self.subTest(filename=filename, expression=expression):
                    self.compare(expression, molecule)


if __name__ == '__main__':
    unittest.main()
