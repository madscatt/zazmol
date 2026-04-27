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

from unittest import main, skipIf
import unittest
import sasmol.system as system
import sasmol.operate as operate

import numpy

import warnings

import os
floattype = os.environ['SASMOL_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'pdb_common') + os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','operate')+os.path.sep

import sasmol.system as system

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'data', 'pdb_common') + os.path.sep

class Test_Operate_Move_Align_PMI(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings('ignore')
        self.o1 = system.Molecule(0)
        self.o1.read_pdb(DataPath + '1CRN.pdb')

    def assert_list_almost_equal(self, a, b, places=5):
        if len(a) != len(b):
            raise TypeError
        else:
            for i in range(len(a)):
                if isinstance(a[i], (int, float, numpy.generic)):
                    if numpy.isnan(a[i]) and numpy.isnan(b[i]):
                        continue
                    self.assertAlmostEqual(a[i], b[i], places)
                else:
                    self.assert_list_almost_equal(a[i], b[i], places)

    def assert_axis_aligned(self, result_axis, expected_axis, places=2):
        result_axis = result_axis / numpy.linalg.norm(result_axis)
        expected_axis = expected_axis / numpy.linalg.norm(expected_axis)
        alignment = abs(numpy.dot(result_axis, expected_axis))
        self.assertAlmostEqual(alignment, 1.0, places)

    def test_align_pmi_on_axis_x(self):
        frame = 0
        pmi_eigenvector = 0
        alignment_vector_axis = 'x'
        self.o1.align_pmi_on_axis(frame, pmi_eigenvector, alignment_vector_axis)
        expected_axis = numpy.array([1.0, 0.0, 0.0])
        uk, ak, I = self.o1.calculate_principal_moments_of_inertia(frame)
        result_axis = ak[:, pmi_eigenvector]
        self.assert_axis_aligned(result_axis, expected_axis, places=2)

    def test_align_pmi_on_axis_y(self):
        frame = 0
        pmi_eigenvector = 1
        alignment_vector_axis = 'y'
        self.o1.align_pmi_on_axis(frame, pmi_eigenvector, alignment_vector_axis)
        expected_axis = numpy.array([0.0, 1.0, 0.0])
        uk, ak, I = self.o1.calculate_principal_moments_of_inertia(frame)
        result_axis = ak[:, pmi_eigenvector]
        self.assert_axis_aligned(result_axis, expected_axis, places=2)

    def test_align_pmi_on_axis_z(self):
        frame = 0
        pmi_eigenvector = 2
        alignment_vector_axis = 'z'
        self.o1.align_pmi_on_axis(frame, pmi_eigenvector, alignment_vector_axis)
        expected_axis = numpy.array([0.0, 0.0, 1.0])
        uk, ak, I = self.o1.calculate_principal_moments_of_inertia(frame)
        result_axis = ak[:, pmi_eigenvector]
        self.assert_axis_aligned(result_axis, expected_axis, places=2)

    def tearDown(self):
        pass

if __name__ == '__main__':
    main()
