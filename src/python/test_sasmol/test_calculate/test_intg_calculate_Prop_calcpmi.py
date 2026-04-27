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
import unittest
import sasmol.system as system

import numpy

import os
floattype = os.environ['SASMOL_FLOATTYPE']

PdbPath = os.path.join(os.path.dirname(os.path.realpath(
    __file__)), '..', 'data', 'pdb_common')+os.path.sep


class Test_sascalc_Prop_calcpmi(unittest.TestCase):

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

    def assert_list_almost_equal_flip_sign_allowed(self, a, b, places=5):
        if (len(a) != len(b)):
            raise TypeError
        else:
            sign = 1
            for i in range(len(a)):
                if isinstance(a[i], (int, float)):
                    if (numpy.isnan(a[i]) and numpy.isnan(b[i])):
                        continue
                    if (a[i]*b[i] < 0.0):
                        sign = -1
                    self.assertAlmostEqual(a[i], sign*b[i], places)
                else:
                    self.assert_list_almost_equal_flip_sign_allowed(
                        a[i], b[i], places)

    def reorder_eigens(self, result_eigenvalues, result_eigenvectors):
        idx = result_eigenvalues.argsort()
        idx = idx[::-1]
        result_eigenvalues = result_eigenvalues[idx]
        result_eigenvectors = result_eigenvectors[idx]
        result_eigenvectors[2] *= -1
        return result_eigenvalues, result_eigenvectors

    def test_null(self):
        with self.assertRaises(Exception):
            self.o.read_pdb(PdbPath+'NULL.pdb')
        with self.assertRaises(Exception):
            result_pmi = self.o.calculate_principal_moments_of_inertia(0)

    def test_one_atom_pdb(self):
        return
        '''
         
        self.o.read_pdb(PdbPath+'1ATM.pdb')
        self.o.calculate_mass()
        result = self.o.calculate_principal_moments_of_inertia(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        expected_I = numpy.zeros((3,3),floattype)
        expected_eigenvalues = numpy.zeros(3,floattype)
        expected_eigenvectors = numpy.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 5)        
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)

        '''

    def test_two_aa_pdb(self):
        self.o.read_pdb(PdbPath+'2AAD.pdb')
        self.o.calculate_mass()
        result = self.o.calculate_principal_moments_of_inertia(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(
            result_eigenvalues, result_eigenvectors)
        # Expected values reflect float32 coordinate storage with float64
        # calculation intermediates.
        expected_I = numpy.array([numpy.array([589.533746311651,  -64.328461565902, -439.375385704833]), numpy.array(
            [-64.328461565902, 1532.135608483722,  -65.398994295632]), numpy.array([-439.375385704833,   -65.398994295632, 1407.583009456943])], floattype)
        expected_eigenvalues = numpy.array(
            [1614.281458831898, 1523.060334881807,  391.910570545354], floattype)
        # expected_eigenvectors = numpy.array([numpy.array([ 0.33751536,  0.41020629, -0.84723915]), numpy.array([ 0.22717091, -0.90894656, -0.3495848 ]), numpy.array([-0.913497  , -0.07447786, -0.39997036])],floattype)
        expected_eigenvectors = numpy.array([numpy.array([0.33751523,  0.41020703, -0.84723885]), numpy.array(
            [0.22717101, -0.90894624, -0.34958556]), numpy.array([-0.91349702, -0.07447765, -0.39997034])], floattype)
        # print("I = ", result_I)
        # print("result eigenvalues  = ", result_eigenvalues)
        # print("result eigenvectors  = ", result_eigenvectors)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_I, result_I, 5)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvalues, result_eigenvalues, 3)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvectors, result_eigenvectors, 3)

    def test_rna_pdb(self):
        self.o.read_pdb(PdbPath+'rna.pdb')
        self.o.calculate_mass()
        result = self.o.calculate_principal_moments_of_inertia(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(
            result_eigenvalues, result_eigenvectors)
        # Expected values reflect float32 coordinate storage with float64
        # calculation intermediates.
        expected_I = numpy.array([numpy.array([3.044118942442e+08,   4.043337055251e+07,   4.065207099889e+07]), numpy.array(
            [4.043337055251e+07,   3.085291059812e+08,  -4.593367594306e+07]), numpy.array([4.065207099889e+07,  -4.593367594306e+07,   3.021965836672e+08])], floattype)
        expected_eigenvalues = numpy.array(
            [3.516875335607e+08, 3.431749508626e+08, 2.202750994692e+08], floattype)
        expected_eigenvectors = numpy.array([numpy.array([-0.1525973, -0.78373478,  0.60205802]), numpy.array(
            [0.8122177,  0.24761209,  0.52819567]), numpy.array([0.56304216, -0.56960341, -0.59877832])], floattype)
        # print(f"{numpy.array2string(result_eigenvalues, precision=12, floatmode='fixed')}")
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_I, result_I, -1)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvalues, result_eigenvalues, 3)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvectors, result_eigenvectors, 3)

    def test_1CRN_pdb(self):
        self.o.read_pdb(PdbPath+'1CRN.pdb')
        self.o.calculate_mass()
        result = self.o.calculate_principal_moments_of_inertia(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(
            result_eigenvalues, result_eigenvectors)
        # Expected values reflect float32 coordinate storage with float64
        # calculation intermediates.
        expected_I = numpy.array([numpy.array([258450.066403634846, -45258.199576749066, -67627.065565671233]), numpy.array([-45258.199576749066,
                                 311061.104958547745,  -3089.063300746677]), numpy.array([-67627.065565671233,  -3089.063300746676, 244380.695794588071])], floattype)
        expected_eigenvalues = numpy.array(
            [349987.999953595398, 288718.592304701044, 175185.274898475007], floattype)
        expected_eigenvectors = numpy.array([numpy.array([0.61882991, -0.68963421, -0.37610398]), numpy.array(
            [-0.37918551, -0.68156978,  0.62584422]), numpy.array([-0.68794469, -0.24467794, -0.68327506])], floattype)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_I, result_I, 3)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvalues, result_eigenvalues, 3)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvectors, result_eigenvectors, 3)

    @skipIf(os.environ['SASMOL_LARGETEST'] == 'n', "I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(PdbPath+'1KP8.pdb')
        self.o.calculate_mass()
        result = self.o.calculate_principal_moments_of_inertia(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(
            result_eigenvalues, result_eigenvectors)
        # Expected values reflect float32 coordinate storage with float64
        # calculation intermediates.
        expected_I = numpy.array([numpy.array([2.118857176574e+09, -5.073118355158e+06, -6.581597703808e+06]), numpy.array(
            [-5.073118355158e+06, 2.118487346086e+09, 7.279001798601e+06]), numpy.array([-6.581597703808e+06, 7.279001798601e+06, 1.903425054207e+09])], floattype)
        expected_eigenvalues = numpy.array(
            [2.124183019558e+09, 2.113597827755e+09, 1.902988729553e+09], floattype)
        expected_eigenvectors = numpy.array([numpy.array([0.717208918858, -0.695447830514, -0.044313448795]), numpy.array(
            [-0.696225741993, -0.717816386552, -0.003058003247]), numpy.array([-0.029682237286, 0.033045390311, -0.999012996396])], floattype)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_I, result_I, -2)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvalues, result_eigenvalues, 3)
        self.assert_list_almost_equal_flip_sign_allowed(
            expected_eigenvectors, result_eigenvectors, 3)

    def test_problem_pdb(self):
        with self.assertRaises(Exception):
            self.o.read_pdb(PdbPath+'1PSI.pdb')
        with self.assertRaises(Exception):
            result_pmi = self.o.calculate_principal_moments_of_inertia(0)

    def tearDown(self):
        pass


if __name__ == '__main__':
    main()
