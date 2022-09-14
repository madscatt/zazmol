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

from sasmol.test_sasmol.utilities import env,util

from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.system as system
import sasmol.operate as operate
import sasmol.calculate as calculate

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASMOL_FLOATTYPE']

class Test_unit_operate_Move_translate(MockerTestCase): 

    def setUp(self):
        self.back_mass_check = operate.Move.mass_check
        self.back_calccom = calculate.Calculate.calculate_center_of_mass 

        self.m = Mocker()

        calculate.Calculate.calculate_center_of_mass = self.m.mock()
        calculate.Calculate.calculate_center_of_mass(ARGS)
        self.m.result(None)
        self.m.count(0,None)

        operate.Move.mass_check = self.m.mock()
        operate.Move.mass_check()
        self.m.result(None)
        self.m.count(0,None)

        self.m.replay()

        self.o=system.Molecule(0)

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if isinstance(a[i],(int,float,numpy.generic)):
                 if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
                 self.assertAlmostEqual(a[i],b[i],places)
              else:
                 self.assert_list_almost_equal(a[i],b[i],places)

    def test_one_atom(self):
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.0]]],floattype))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[0.0, 5.0, 9.0]]], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_two_atoms(self):
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.87],[-5.0, 3.2, 6.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-4.0, 2.0, 9.97],[-8.0, 3.2, 12.1]]], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-1.8, 2.0, 9.1],[-5.0, 5.0, 12.1],[4.0, 8.0, 15.1],[-2.0, 3.0, 11.1],[-1.0, 4.0, 12.1],[-3.0, 2.0, 9.1]]],floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_inf1(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, util.HUGE],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-1.8, 2.0, util.HUGE],[-5.0, 5.0, 12.1],[4.0, 8.0, 15.1],[-2.0, 3.0, 11.1],[-1.0, 4.0, 12.1],[-3.0, 2.0, 9.1]]],floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_inf2(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, util.INF],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-1.8, 2.0, util.INF],[-5.0, 5.0, 12.1],[4.0, 8.0, 15.1],[-2.0, 3.0, 11.1],[-1.0, 4.0, 12.1],[-3.0, 2.0, 9.1]]],floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_nan(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, util.NAN],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-1.8, 2.0, util.NAN],[-5.0, 5.0, 12.1],[4.0, 8.0, 15.1],[-2.0, 3.0, 11.1],[-1.0, 4.0, 12.1],[-3.0, 2.0, 9.1]]],floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_tiny(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, util.TINY],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-1.8, 2.0, util.TINY+6.1],[-5.0, 5.0, 12.1],[4.0, 8.0, 15.1],[-2.0, 3.0, 11.1],[-1.0, 4.0, 12.1],[-3.0, 2.0, 9.1]]],floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_zero(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, util.ZERO],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        value = numpy.array([-3.0, 0.0, 6.1],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print(result_coor)
        expected_coor = numpy.array([[[-1.8, 2.0, 6.1],[-5.0, 5.0, 12.1],[4.0, 8.0, 15.1],[-2.0, 3.0, 11.1],[-1.0, 4.0, 12.1],[-3.0, 2.0, 9.1]]],floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def tearDown(self):
        self.m.verify()
        calculate.Calculate.calculate_center_of_mass=self.back_calccom
        operate.Move.mass_check  =self.back_mass_check


if __name__ == '__main__': 
   main() 

