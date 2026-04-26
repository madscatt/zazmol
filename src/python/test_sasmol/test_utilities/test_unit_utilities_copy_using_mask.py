'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

from sasmol.test_sasmol.utilities import env, util

from unittest import main
import unittest
import sasmol.system as system
import sasmol.utilities as utilities

import numpy
import warnings
import os

floattype = os.environ['SASMOL_FLOATTYPE']


class Test_unit_utilities_copy_using_mask(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings('ignore')
        self.o = system.Molecule(0)

        # realistic minimal molecule
        self.o.setName(['N', 'CA', 'CB', 'CG'])
        self.o.setElement(['C', 'C', 'C', 'C'])
        self.o.setResid([1,1,1,1])
        self.o.setIndex([1,2,3,4])
        self.o.setOriginal_index(numpy.array([1,2,3,4]))
        self.o.setOriginal_resid(numpy.array([1,1,1,1]))

        self.o.setCoor(numpy.array([[
            [1.0,2.0,3.0],
            [4.0,5.0,6.0],
            [7.0,8.0,9.0],
            [10.0,11.0,12.0]
        ]], floattype))

        self.o.setMass(numpy.array([12.0,12.0,12.0,12.0]))
        self.o.setNatoms(4)

    def test_mask_subset(self):
        mask = [0, 2]
        new = utilities.Copy_Using_Mask.from_sasmol(self.o, mask)

        self.assertEqual(len(new.name()), 2)
        self.assertEqual(new.name()[0], 'N')
        self.assertEqual(new.name()[1], 'CB')

    def test_coordinates_subset(self):
        mask = [1, 3]
        new = utilities.Copy_Using_Mask.from_sasmol(self.o, mask)

        coor = new.coor()[0]
        self.assertEqual(coor.shape[0], 2)
        self.assertTrue((coor[0] == numpy.array([4.0,5.0,6.0])).all())
        self.assertTrue((coor[1] == numpy.array([10.0,11.0,12.0])).all())

    def test_index_consistency(self):
        mask = [0,1,2]
        new = utilities.Copy_Using_Mask.from_sasmol(self.o, mask)

        self.assertEqual(len(new.index()), 3)
        self.assertEqual(new.index()[0], 1)
        self.assertEqual(new.index()[2], 3)

    def tearDown(self):
        pass


if __name__ == '__main__':
    main()
