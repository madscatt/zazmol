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

contract

make sure the keys are unique
make sure the right protein list was generated
'''

from unittest import main
import unittest

import warnings

import sasmol.system as system

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', 'data', 'sasmol', 'properties') + os.path.sep


class Test_unit_properties_Atomic_protein_scattering_lengths(
        unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings('ignore')
        self.o = system.Molecule(0)
        self.protein_scattering_lengths = (
            self.o.protein_scattering_lengths())

    def unique(self, seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if x not in seen and not seen_add(x)]

    def test_uniqe(self):
        protein_scattering_lengths_unique = self.unique(
            list(self.protein_scattering_lengths.keys()))
        self.assertEqual(
            list(self.protein_scattering_lengths.keys()),
            protein_scattering_lengths_unique)

    def test_all(self):
        datafile = DataPath + 'protein_scattering_lengths.txt'
        ele = {}
        for line in open(datafile).readlines():
            aa = line.split()[0]
            ele[aa] = eval(line[4:])
        self.assertEqual(ele, self.protein_scattering_lengths)


if __name__ == '__main__':
    main()
