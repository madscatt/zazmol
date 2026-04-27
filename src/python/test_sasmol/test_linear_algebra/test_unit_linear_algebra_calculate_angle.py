from sasmol.test_sasmol.utilities import env, util

from unittest import main
import unittest
import sasmol.linear_algebra as linear_algebra
import numpy
import warnings
import math

import os
floattype = os.environ['SASMOL_FLOATTYPE']


class Test_linear_algebra_calculate_angle(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings('ignore')

    def test_zero_vectors(self):
        a = numpy.array([0.0, 0.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([0.0, 0.0, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_orthogonal(self):
        a = numpy.array([1.0, 0.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([0.0, 1.0, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        expected = math.pi / 2.0
        self.assertAlmostEqual(result, expected)

    def test_parallel(self):
        a = numpy.array([1.0, 0.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([2.0, 0.0, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        expected = 0.0
        self.assertAlmostEqual(result, expected)

    def test_antiparallel(self):
        a = numpy.array([1.0, 0.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([-1.0, 0.0, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        expected = math.pi
        self.assertAlmostEqual(result, expected)

    def test_arbitrary(self):
        a = numpy.array([1.0, 2.0, 3.0], floattype)
        b = numpy.array([-1.0, 6.0, 8.0], floattype)
        c = numpy.array([-4.0, -1.0, 4.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)

        # expected from direct formula
        u = a - b
        v = c - b
        expected = math.acos(numpy.dot(u, v) / (numpy.linalg.norm(u) * numpy.linalg.norm(v)))

        self.assertAlmostEqual(result, expected)

    def test_inf(self):
        a = numpy.array([util.INF, 1.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([1.0, 0.0, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_nan(self):
        a = numpy.array([util.NAN, 1.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([1.0, 0.0, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        self.assertTrue(numpy.isnan(result))

    def test_tiny(self):
        a = numpy.array([util.TINY, 1.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([1.0, util.TINY, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        expected = math.pi / 2.0
        self.assertAlmostEqual(result, expected, places=5)

    def test_zero(self):
        a = numpy.array([util.ZERO, 1.0, 0.0], floattype)
        b = numpy.array([0.0, 0.0, 0.0], floattype)
        c = numpy.array([1.0, util.ZERO, 0.0], floattype)
        result = linear_algebra.calculate_angle(a, b, c)
        expected = math.pi / 2.0
        self.assertAlmostEqual(result, expected, places=5)

    def tearDown(self):
        pass


if __name__ == '__main__':
    main()
