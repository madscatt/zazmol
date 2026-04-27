'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
from sasmol.test_sasmol.utilities import env

from unittest import main, skipIf
import sasmol.system as system
import unittest

import numpy

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', 'data', 'pdb_common')+os.path.sep
dcdDataPath = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', 'data', 'dcd_common')+os.path.sep


def hugeDcdPath(filename):
    return os.path.join('/tmp', filename)


class Test_sascalc_Prop_calc_minmax_all_steps(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def assert_list_almost_equal(self, a, b, places=5):
        if (len(a) != len(b)):
            raise TypeError
        else:
            for i in range(len(a)):
                if (numpy.isnan(a[i]) and numpy.isnan(b[i])):
                    continue
                self.assertAlmostEqual(a[i], b[i], places)

    def test_null(self):
        with self.assertRaises(Exception):
            self.o.calcminmax_frame(0)

    def test_one_atom(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            dcdDataPath+'1ATM.dcd')
        expected_minmax = [[73.944, 38.799, 41.652], [76.944,  41.799, 41.652]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_two_aa(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            dcdDataPath+'2AAD.dcd')
        expected_minmax = [[-79.712, -46.273,  39.354],
                           [79.712,  46.273,  43.910]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_rna_1to10(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            dcdDataPath+'rna-1to10.dcd')
        expected_minmax = [[-43.801, -44.888, -42.605],
                           [41.234,  39.706,  41.903]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_calc_minmax_all_steps_alias_matches_primary_method(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')

        primary = self.o.calculate_minimum_and_maximum_all_steps(
            dcdDataPath+'1ATM.dcd')
        alias = self.o.calc_minmax_all_steps(dcdDataPath+'1ATM.dcd')

        self.assert_list_almost_equal(primary[0], alias[0])
        self.assert_list_almost_equal(primary[1], alias[1])

    @skipIf(os.environ['SASMOL_HUGETEST'] == 'n', "I am not testing huge files")
    def test_rna_0point8g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            hugeDcdPath('rna-0.8g.dcd'))
        expected_minmax = [[-88.148, -86.246, -81.494],
                           [84.491, 77.158, 84.429]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)

    @skipIf(os.environ['SASMOL_HUGETEST'] == 'n', "I am not testing huge files")
    def test_rna_1point0g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            hugeDcdPath('rna-1.0g.dcd'))
        expected_minmax = [[-88.148, -86.246, -81.494],
                           [84.491, 77.158, 84.429]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)

    @skipIf(os.environ['SASMOL_HUGETEST'] == 'n', "I am not testing huge files")
    def test_rna_2point0g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            hugeDcdPath('rna-2.0g.dcd'))
        expected_minmax = [[-88.148, -86.246, -81.494],
                           [84.491, 77.158, 84.429]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)

    @skipIf(os.environ['SASMOL_HUGETEST'] == 'n', "I am not testing huge files")
    def test_rna_3point2g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            hugeDcdPath('rna-3.2g.dcd'))
        expected_minmax = [[-88.148, -86.246, -81.494],
                           [84.491, 77.158, 84.429]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)

    @skipIf(os.environ['SASMOL_HUGETEST'] == 'n', "I am not testing huge files")
    def test_rna_6point4g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax = self.o.calculate_minimum_and_maximum_all_steps(
            hugeDcdPath('rna-6.4g.dcd'))
        expected_minmax = [[-88.148, -86.246, -81.494],
                           [84.491, 77.158, 84.429]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1], 3)

    def tearDown(self):
        pass


if __name__ == '__main__':
    main()
