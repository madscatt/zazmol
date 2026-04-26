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

import unittest
from unittest import mock

import numpy

import sasmol.system as system
import sasmol.view_vmd as view_vmd


class Test_intg_system_Atom_send_coordinates_to_vmd(unittest.TestCase):

    @mock.patch('sasmol.view.view_vmd.send_coordinates_to_vmd')
    def test_send_coordinates_to_vmd_passes_float32_coordinates(self, mock_send):
        molecule = system.Molecule(0)
        molecule.setCoor(numpy.array(
            [[[1.25, 2.5, 3.75], [4.0, 5.5, 6.75]]], numpy.float64
        ))

        molecule.send_coordinates_to_vmd(1085, False)

        mock_send.assert_called_once()
        tx, ty, tz, port, flag = mock_send.call_args[0]

        self.assertEqual(port, 1085)
        self.assertFalse(flag)
        self.assertEqual(tx.dtype, numpy.float32)
        self.assertEqual(ty.dtype, numpy.float32)
        self.assertEqual(tz.dtype, numpy.float32)
        self.assertTrue(numpy.array_equal(
            tx, numpy.array([1.25, 4.0], numpy.float32)))
        self.assertTrue(numpy.array_equal(
            ty, numpy.array([2.5, 5.5], numpy.float32)))
        self.assertTrue(numpy.array_equal(
            tz, numpy.array([3.75, 6.75], numpy.float32)))

    def test_view_vmd_rejects_mismatched_array_lengths(self):
        x = numpy.array([1.0], numpy.float32)
        y = numpy.array([2.0, 3.0], numpy.float32)
        z = numpy.array([4.0], numpy.float32)

        with self.assertRaisesRegex(ValueError, r"Arrays of lengths \(1,2,1\) given"):
            view_vmd.send_coordinates_to_vmd(x, y, z, 1085, False)

    def test_view_vmd_rejects_non_float32_arrays(self):
        x = numpy.array([1.0], numpy.float64)
        y = numpy.array([2.0], numpy.float64)
        z = numpy.array([3.0], numpy.float64)

        with self.assertRaisesRegex(TypeError, r"Coordinate arrays must have dtype float32"):
            view_vmd.send_coordinates_to_vmd(x, y, z, 1085, False)


if __name__ == '__main__':
    unittest.main()
