import unittest

import numpy

import sasmol.config as config
import sasmol.system as system
import sasmol.utilities as utilities


class Test_subset_Mask_duplicate_molecule(unittest.TestCase):

    def make_molecule(self):
        molecule = system.Molecule(0)
        molecule.setNatoms(2)
        molecule.setCoor(numpy.arange(12, dtype=config.COORD_DTYPE).reshape(2, 2, 3))
        molecule.setName(['CA', 'N'])
        molecule.setResid(numpy.array([1, 2], dtype=numpy.int32))
        molecule.setResname(['ALA', 'GLY'])
        molecule.setChain(['A', 'B'])
        molecule.setSegname(['SEG1', 'SEG2'])
        molecule.setHeader(['REMARK test metadata'])
        molecule.setConect([[1, 2]])
        molecule.setBiomt({1: {
            'subdivs': ['A', 'B'],
            'rot': [numpy.eye(3)],
            'trans': [numpy.zeros(3)],
        }})
        return molecule

    def assert_initial_values(self, molecule):
        expected = self.make_molecule()
        self.assertEqual(molecule.natoms(), 2)
        self.assertEqual(molecule.coor().shape, (2, 2, 3))
        self.assertEqual(molecule.coor().dtype, numpy.dtype(config.COORD_DTYPE))
        numpy.testing.assert_array_equal(molecule.coor(), expected.coor())
        numpy.testing.assert_array_equal(molecule.resid(), expected.resid())
        self.assertEqual(molecule.resid().dtype, expected.resid().dtype)
        for field in ('name', 'resname', 'chain', 'segname', 'header', 'conect'):
            self.assertEqual(getattr(molecule, field)(), getattr(expected, field)())
        self.assertEqual(set(molecule.biomt()), {1})
        metadata = molecule.biomt()[1]
        self.assertEqual(metadata['subdivs'], ['A', 'B'])
        self.assertEqual(len(metadata['rot']), 1)
        self.assertEqual(len(metadata['trans']), 1)
        numpy.testing.assert_array_equal(metadata['rot'][0], numpy.eye(3))
        numpy.testing.assert_array_equal(metadata['trans'][0], numpy.zeros(3))

    def test_empty_molecule_method_returns_one_copy(self):
        molecule = system.Molecule(0)
        copies = molecule.duplicate_molecule(1)
        self.assertIsInstance(copies, list)
        self.assertEqual(len(copies), 1)
        self.assertIsNot(copies[0], molecule)

    def test_both_entry_points_preserve_values_and_deep_copy_independence(self):
        for entry_point in ('method', 'utility'):
            with self.subTest(entry_point=entry_point):
                original = self.make_molecule()
                copies = (original.duplicate_molecule(2) if entry_point == 'method'
                          else utilities.duplicate_molecule(original, 2))
                self.assertIsInstance(copies, list)
                self.assertEqual(len(copies), 2)
                first, sibling = copies
                self.assertIsNot(first, original)
                self.assertIsNot(sibling, original)
                self.assertIsNot(first, sibling)
                for molecule in (original, first, sibling):
                    self.assertIs(type(molecule), system.Molecule)
                    self.assert_initial_values(molecule)

                first.coor()[1, 0, 0] = -99
                first.name()[0] = 'CHANGED'
                first.resid()[0] = 99
                first.segname()[0] = 'CHANGED'
                first.header().append('CHANGED')
                first.conect()[0].append(3)
                first.biomt()[1]['subdivs'].append('C')
                first.biomt()[1]['rot'][0][0, 0] = 5
                first.biomt()[1]['trans'][0][0] = 10
                self.assert_initial_values(original)
                self.assert_initial_values(sibling)

    def test_method_preserves_utility_count_behavior(self):
        molecule = system.Molecule(0)
        for count, expected in ((0, 0), (-1, 0), (False, 0), (True, 1), (2, 2)):
            with self.subTest(count=count):
                copies = molecule.duplicate_molecule(count)
                self.assertIsInstance(copies, list)
                self.assertEqual(len(copies), expected)
                self.assertEqual(len(copies), len(utilities.duplicate_molecule(molecule, count)))
        with self.assertRaises(TypeError):
            molecule.duplicate_molecule(1.5)

    def test_reserved_kwargs_remain_ignored(self):
        self.assertEqual(len(system.Molecule(0).duplicate_molecule(1, future=True)), 1)


if __name__ == '__main__':
    unittest.main()
