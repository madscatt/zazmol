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

import numpy

import sasmol.config as config
import sasmol.system as system


class Test_sascalc_Prop_calcresidue_charge(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def test_two_residues(self):
        self.o.setNatoms(4)
        self.o.setResid([1, 1, 2, 2])
        self.o.setAtom_charge(numpy.array([0.1, 0.2, -0.4, 0.3],
                                          config.CALC_DTYPE))

        self.o.calculate_residue_charge()

        expected = numpy.array([0.3, 0.3, -0.1, -0.1], config.CALC_DTYPE)
        self.assertEqual(self.o.residue_charge().dtype, config.CALC_DTYPE)
        numpy.testing.assert_allclose(self.o.residue_charge(), expected)


if __name__ == '__main__':
    main()
