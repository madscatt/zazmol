'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

import unittest

import numpy

import sasmol.config as config
import sasmol.system as system


class Test_subset_Mask_set_coor_using_mask(unittest.TestCase):

    def make_target(self):
        coor = numpy.array([[
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0],
            [4.0, 4.0, 4.0],
        ]], dtype=config.COORD_DTYPE)
        return system.Molecule_Maker(4, coor=coor)

    def make_source(self):
        coor = numpy.array([[
            [10.0, 10.0, 10.0],
            [20.0, 20.0, 20.0],
        ]], dtype=config.COORD_DTYPE)
        return system.Molecule_Maker(2, coor=coor)

    def test_selected_atoms_are_replaced(self):
        target = self.make_target()
        source = self.make_source()
        mask = numpy.array([0, 1, 0, 1], dtype=numpy.int32)

        error = target.set_coor_using_mask(source, 0, mask)

        self.assertEqual(error, [])
        self.assertTrue(numpy.allclose(target.coor()[0, 1], [10.0, 10.0, 10.0]))
        self.assertTrue(numpy.allclose(target.coor()[0, 3], [20.0, 20.0, 20.0]))

    def test_unselected_atoms_remain_unchanged(self):
        target = self.make_target()
        original = target.coor().copy()
        source = self.make_source()
        mask = numpy.array([0, 1, 0, 1], dtype=numpy.int32)

        error = target.set_coor_using_mask(source, 0, mask)

        self.assertEqual(error, [])
        self.assertTrue(numpy.allclose(target.coor()[0, 0], original[0, 0]))
        self.assertTrue(numpy.allclose(target.coor()[0, 2], original[0, 2]))

    def test_incompatible_source_frame_returns_error(self):
        target = self.make_target()
        source = self.make_source()
        mask = numpy.array([0, 1, 0, 1], dtype=numpy.int32)

        error = target.set_coor_using_mask(source, 1, mask)

        self.assertTrue(len(error) > 0)


if __name__ == '__main__':
    unittest.main()
