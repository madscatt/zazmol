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

if __name__ == '__main__':
    main()
