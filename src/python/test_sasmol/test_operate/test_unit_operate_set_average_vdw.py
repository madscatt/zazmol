'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

from sasmol.test_sasmol.utilities import env, util

from unittest import main
import unittest
import sasmol.system as system

import numpy
import warnings
import os

floattype = os.environ['SASMOL_FLOATTYPE']


class Test_unit_operate_set_average_vdw(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings('ignore')

        self.o = system.Molecule(0)

        # minimal realistic molecule
        self.o.setElement(['C', 'N', 'O'])
        self.o.setNatoms(3)

        self.o.setCoor(numpy.array([[
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0]
        ]], floattype))

    def test_set_average_vdw_basic(self):
        self.o.set_average_vdw()

        vdw = self.o.atom_vdw()

        # length check
        self.assertEqual(len(vdw), 3)

        # structure check (natoms, 2)
        for entry in vdw:
            self.assertEqual(len(entry), 2)

        # radii are in second column
        self.assertAlmostEqual(vdw[0][1], 2.00, places=2)
        self.assertAlmostEqual(vdw[1][1], 1.85, places=2)
        self.assertAlmostEqual(vdw[2][1], 1.74, places=2)

    def test_unknown_element(self):
        self.o.setElement(['X', 'C', 'N'])
        self.o.setNatoms(3)

        self.o.set_average_vdw()
        vdw = self.o.atom_vdw()

        # unknown element handled but structure preserved
        self.assertEqual(len(vdw[0]), 2)
        self.assertEqual(len(vdw[1]), 2)

    def test_length_consistency(self):
        self.o.set_average_vdw()
        self.assertEqual(len(self.o.atom_vdw()), self.o.natoms())

    def test_structure(self):
        self.o.set_average_vdw()
        vdw = self.o.atom_vdw()

        for entry in vdw:
            self.assertEqual(len(entry), 2)


    def tearDown(self):
        pass


if __name__ == '__main__':
    main()
