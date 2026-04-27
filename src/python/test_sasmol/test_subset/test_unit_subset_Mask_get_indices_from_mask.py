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

from unittest import main
import unittest

import sasmol.system as system


class Test_unit_subset_Mask_get_indices_from_mask(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def run_case(self, natoms, mask, expected_indices):
        self.o.setNatoms(natoms)
        result_indices = self.o.get_indices_from_mask(mask)
        self.assertEqual(list(result_indices), expected_indices)

    def test_null(self):
        self.run_case(0, [], [])

    def test_1atom_mask_none(self):
        self.run_case(1, [0], [])

    def test_1atom_mask_one(self):
        self.run_case(1, [1], [0])

    def test_2atom_mask_none(self):
        self.run_case(2, [0, 0], [])

    def test_2atom_mask_first(self):
        self.run_case(2, [1, 0], [0])

    def test_2atom_mask_second(self):
        self.run_case(2, [0, 1], [1])

    def test_2atom_mask_both(self):
        self.run_case(2, [1, 1], [0, 1])

    def test_1000atoms_mask_none(self):
        self.run_case(1000, [0] * 1000, [])

    def test_1000atoms_mask_even(self):
        self.run_case(1000, [1, 0] * 500, list(range(0, 1000, 2)))

    def test_1000atoms_mask_all(self):
        self.run_case(1000, [1] * 1000, list(range(1000)))

    def test_negative_longer_mask(self):
        self.o.setNatoms(10)
        with self.assertRaises(Exception):
            self.o.get_indices_from_mask([1] * 20)

    def test_negative_shorter_mask(self):
        self.o.setNatoms(10)
        with self.assertRaises(Exception):
            self.o.get_indices_from_mask([1] * 5)

    def test_scalar_mask_selects_all_atoms(self):
        self.o.setNatoms(10)
        self.assertEqual(
            list(self.o.get_indices_from_mask(1)),
            list(range(10)))


if __name__ == '__main__':
    main()
