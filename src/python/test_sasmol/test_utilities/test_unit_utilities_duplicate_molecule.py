from unittest import main
import unittest
import sasmol.system as system
import sasmol.utilities as utilities

class Test(unittest.TestCase):
    def test_dup(self):
        o = system.Molecule(0)
        mols = utilities.duplicate_molecule(o, 2)
        self.assertEqual(len(mols), 2)
        self.assertNotEqual(id(mols[0]), id(mols[1]))
        self.assertIsNot(mols[0], o)
        self.assertIsNot(mols[1], o)

    def test_count_behavior(self):
        molecule = system.Molecule(0)
        for count, expected in ((0, 0), (-1, 0), (False, 0), (True, 1), (2, 2)):
            with self.subTest(count=count):
                copies = utilities.duplicate_molecule(molecule, count)
                self.assertIsInstance(copies, list)
                self.assertEqual(len(copies), expected)
        with self.assertRaises(TypeError):
            utilities.duplicate_molecule(molecule, 1.5)

if __name__ == '__main__':
    main()
