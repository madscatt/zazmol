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

from sasmol.test_sasmol.utilities import env, util

from unittest import main
import unittest
from sasmol.linear_algebra import matrix_multiply  # Import the C extension method

import numpy

import os
floattype = os.environ.get('SASMOL_FLOATTYPE', 'float32')
#floattype = numpy.float32

class TestMatrixMultiply(unittest.TestCase):

    def assert_list_almost_equal(self, a, b, places=6):
        for i in range(len(a)):
            if numpy.isnan(a[i]) and numpy.isnan(b[i]):
                continue    
            self.assertAlmostEqual(a[i], b[i], places=places)

    def test_all_zero_arrays_1(self):
        a = numpy.zeros((3, 3), dtype=floattype)
        b = numpy.zeros((3, 3), dtype=floattype)
        result_error, result = matrix_multiply(a, b)
        expected = numpy.zeros((3, 3), dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_all_zero_arrays_2(self):
        a = numpy.zeros((2, 3), dtype=floattype)
        b = numpy.zeros((3, 2), dtype=floattype)
        result_error, result = matrix_multiply(a, b)
        expected = numpy.zeros((2, 2), dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_all_zero_arrays_3(self):
        a = numpy.zeros((3, 2), dtype=floattype)
        b = numpy.zeros((2, 3), dtype=floattype)
        result_error, result = matrix_multiply(a, b)
        expected = numpy.zeros((3, 3), dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_unit_arrays_1(self):
        a = numpy.array([[1.0, 0.0, 0.0]], dtype=floattype)
        b = numpy.array([[1.0, 1.0, 0.0]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[1.0]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_unit_arrays_2(self):
        a = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=floattype)
        b = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[1.0, 0.0], [0.0, 1.0]], dtype=floattype)

        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_unit_arrays_3(self):
        a = numpy.array([[1.0, 0.0, 0.0]], dtype=floattype)
        b = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[1.0, 0.0]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_arb_1(self):
        a = numpy.array([[12.0, -20.0, -80.0], [2.02, -901.0, 0.0]], dtype=floattype)
        b = numpy.array([[1.23, 0.03, 20.0], [10.0, 1.0, 3.0]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[-1585.84, -140.0], [-24.5454, -880.8]], dtype=floattype)
        result = numpy.round(result, 3)
        expected = numpy.round(expected, 3)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_arb_2(self):
        a = numpy.array([[1.0168308, -1.35572028, -1.35362422], [-0.69958848, 1.66901076, 0.49978462]], dtype=floattype)
        b = numpy.array([[-20.69958848, 16.66901076, 20.49978462]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[-71.396], [52.547]], dtype=floattype)
        result = numpy.round(result, 3)
        expected = numpy.round(expected, 3)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_inf_1(self):
        a = numpy.array([[util.INF, -20.0, 80.0]], dtype=floattype)
        b = numpy.array([[1.23, 2.0, 20.0]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[util.INF]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_inf_2(self):
        a = numpy.array([[util.INF, -20.0, 80.0], [2.02, util.INF, 20.0]], dtype=floattype)
        b = numpy.array([[1.23, 2.0, 20.0], [10.0, 21.0, util.INF]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[util.INF, util.INF], [util.INF, util.INF]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_nan(self):
        a = numpy.array([[util.NAN, -20.0, 80.0], [2.02, util.NAN, 20.0]], dtype=floattype)
        b = numpy.array([[1.23, 2.0, 20.0], [10.0, 21.0, util.NAN]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[util.NAN, util.NAN], [util.NAN, util.NAN]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_nan_inf(self):
        a = numpy.array([[util.INF, -20.0, 80.0], [2.02, util.NAN, 20.0]], dtype=floattype)
        b = numpy.array([[1.23, 2.0, 20.0], [10.0, 21.0, util.NAN]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[util.INF, util.NAN], [util.NAN, util.NAN]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_tiny(self):
        a = numpy.array([[util.TINY, -20.0, 80.0], [2.02, util.TINY, 20.0]], dtype=floattype)
        b = numpy.array([[1.23, 2.0, 20.0], [10.0, 21.0, util.TINY]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[1560.0, -420.0], [402.4846, 20.2]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

    def test_zero(self):
        a = numpy.array([[util.ZERO, -20.0, 80.0], [2.02, util.ZERO, 20.0]], dtype=floattype)
        b = numpy.array([[1.23, 2.0, 20.0], [10.0, 21.0, util.ZERO]], dtype=floattype).T
        result_error, result = matrix_multiply(a, b)
        expected = numpy.array([[1560.0, -420.0], [402.4846, 20.2]], dtype=floattype)
        self.assert_list_almost_equal(result.flatten(), expected.flatten())
        expected_error = []
        self.assertEqual(result_error, expected_error)

if __name__ == '__main__':
    unittest.main()
