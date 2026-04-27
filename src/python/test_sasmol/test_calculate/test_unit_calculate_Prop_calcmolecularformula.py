'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

import os
import unittest

import sasmol.system as system


DATA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '..',
    'data',
    'pdb_common'
) + os.path.sep


class Test_sascalc_Prop_calcmolecularformula(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def test_returns_expected_formula_for_known_fixture(self):
        self.o.read_pdb(DATA_PATH + '1ATM.pdb')

        result = self.o.calculate_molecular_formula()

        self.assertEqual(result, {'N': 1})

    def test_updates_formula_descriptor_consistently(self):
        self.o.read_pdb(DATA_PATH + '2AAD.pdb')

        result = self.o.calculate_molecular_formula()

        self.assertEqual(self.o.formula(), result)
        self.assertEqual(result['N'], 2)
        self.assertEqual(result['O'], 3)
        self.assertEqual(result['C'], 10)


if __name__ == '__main__':
    unittest.main()
