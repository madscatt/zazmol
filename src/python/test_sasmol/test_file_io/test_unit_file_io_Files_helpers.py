'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

import unittest

import numpy

import sasmol.config as config
import sasmol.system as system


class Test_unit_file_io_Files_helpers(unittest.TestCase):

    def test_check_for_all_zero_columns_nudges_first_atom_only(self):
        molecule = system.Molecule_Maker(2)
        coor = numpy.array([[
            [0.0, 1.0, 2.0],
            [0.0, 3.0, 4.0],
        ]], dtype=config.COORD_DTYPE)

        molecule.check_for_all_zero_columns(coor)

        self.assertAlmostEqual(float(coor[0, 0, 0]), 1.0e-10, places=15)
        self.assertEqual(float(coor[0, 1, 0]), 0.0)
        self.assertEqual(float(coor[0, 0, 1]), 1.0)
        self.assertEqual(float(coor[0, 1, 2]), 4.0)

    def test_check_for_all_zero_columns_leaves_nonzero_columns_unchanged(self):
        molecule = system.Molecule_Maker(2)
        coor = numpy.array([[
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
        ]], dtype=config.COORD_DTYPE)
        original = coor.copy()

        molecule.check_for_all_zero_columns(coor)

        self.assertTrue(numpy.array_equal(coor, original))

    def test_create_conect_pdb_lines_remaps_and_sorts_indices(self):
        molecule = system.Molecule(0)
        molecule.setOriginal_index(numpy.array([10, 20, 30], numpy.int32))
        molecule.setIndex(numpy.array([1, 5, 3], numpy.int32))
        molecule.setConect({10: [20, 30], 30: [10]})

        conect_lines = molecule.create_conect_pdb_lines()

        self.assertEqual(conect_lines, [
            'CONECT    1    5    3',
            'CONECT    3    1',
        ])


if __name__ == '__main__':
    unittest.main()
