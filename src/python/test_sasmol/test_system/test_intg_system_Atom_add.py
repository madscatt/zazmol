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

import os
import unittest

import numpy

import sasmol.system as system


PDB_DATA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'pdb_common'
) + os.path.sep


class Test_intg_system_Atom_add(unittest.TestCase):

    def test_add_combines_descriptors_and_coordinates(self):
        left = system.Molecule(0)
        right = system.Molecule(1)

        left.read_pdb(PDB_DATA_PATH + '1ATM.pdb')
        right.read_pdb(PDB_DATA_PATH + '1ATM.pdb')

        expected_atom = left.atom() + right.atom()
        expected_name = left.name() + right.name()
        expected_resid = list(left.resid()) + list(right.resid())
        expected_coor = numpy.concatenate((left.coor(), right.coor()), axis=1)

        result = left + right

        self.assertIsNone(result)
        self.assertEqual(left.natoms(), 2)
        self.assertEqual(left.atom(), expected_atom)
        self.assertEqual(left.name(), expected_name)
        self.assertEqual(list(left.index()), [1, 2])
        self.assertEqual(list(left.resid()), expected_resid)
        self.assertTrue(numpy.array_equal(left.coor(), expected_coor))


if __name__ == '__main__':
    unittest.main()
