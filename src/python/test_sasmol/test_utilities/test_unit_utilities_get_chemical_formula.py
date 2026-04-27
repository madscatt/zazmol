'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

from sasmol.test_sasmol.utilities import env, util

import unittest
import sasmol.utilities as utilities

import warnings; warnings.filterwarnings('ignore')


class Test_unit_utilities_get_chemical_formula(unittest.TestCase):

    def test_h2o(self):
        err, formula = utilities.get_chemical_formula("H2O")

        self.assertEqual(err, [])
        self.assertEqual(formula['H'], 2)
        self.assertEqual(formula['O'], 1)

    def test_glucose(self):
        err, formula = utilities.get_chemical_formula("C6H12O6")

        self.assertEqual(err, [])
        self.assertEqual(formula['C'], 6)
        self.assertEqual(formula['H'], 12)
        self.assertEqual(formula['O'], 6)

    def test_parentheses(self):
        err, formula = utilities.get_chemical_formula("Ca(OH)2")

        self.assertEqual(err, [])
        self.assertEqual(formula['Ca'], 1)
        self.assertEqual(formula['O'], 2)
        self.assertEqual(formula['H'], 2)

    def test_nested(self):
        err, formula = utilities.get_chemical_formula("(NH4)2SO4")

        self.assertEqual(err, [])
        self.assertEqual(formula['N'], 2)
        self.assertEqual(formula['H'], 8)
        self.assertEqual(formula['S'], 1)
        self.assertEqual(formula['O'], 4)

    def test_invalid(self):
        err, formula = utilities.get_chemical_formula("XyZ123")

        self.assertTrue(len(err) > 0)


if __name__ == '__main__':
    unittest.main()
