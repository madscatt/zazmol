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

from unittest import main, skipIf
import unittest
import sasmol.system as system
import sasmol.operate as operate

import numpy

import warnings

import os
floattype=os.environ['SASMOL_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_intg_operate_Move_align(unittest.TestCase): 

    def setUp(self):
        warnings.filterwarnings('ignore')
        self.o1=system.Molecule(0)
        self.o2=system.Molecule(1)

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


    #def test_rna_pdb(self):
    def test_1CRN_pdb(self):
        '''

        read in the same PDB file into two separate molecules

        rotate the second molecule and 

        set a basis then align the second molecule on the first molecule

        '''

        print(DataPath+'1CRN.pdb')

        self.o1.read_pdb(DataPath+'1CRN.pdb')
        self.o2.read_pdb(DataPath+'1CRN.pdb')
        
        axis = 'z'
        frame = 0
        theta=numpy.pi/2.0
        
        self.o2.write_pdb("test.pdb", frame, "w")

        self.o2.rotate(frame,axis,theta)

        self.o2.write_pdb("test_rotated.pdb", frame, "w")
      
        basis_1 = 'name[i] == "CA" and (resid[i] >= 20 and resid[i] <= 31)'
        basis_2 = 'name[i] == "CA" and (resid[i] >= 20 and resid[i] <= 31)'

        self.o2.align(self.o1, basis_1, basis_2)

        expected_com = self.o1.calculate_center_of_mass(frame)

        result_com = self.o2.calculate_center_of_mass(frame)

        print('expected_com : calculated = ', expected_com)

        print('result_com : of aligned mol', result_com)

        #
        expected_com = numpy.array([10.21108263, -9.38440245,  6.97775639], floattype)

        self.o2.write_pdb("test_rotated_then_aligned.pdb", frame, "w")

        self.assert_list_almost_equal(expected_com, result_com, 2)

    ### methods below are from rotate test, they have not been altered.

    '''

    def test_1CRN_pdb_rotate(self):
    #def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        axis = 'z'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.rotate(frame,axis,theta)
        result_com  = self.o.calculate_center_of_mass(0)
        #
        expected_com = numpy.array([-9.775, 9.300, 6.978], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)


    @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        axis = 'x'
        frame = 0
        theta=12.0
        #
        self.o.rotate(frame,axis,theta)
        result_com  = self.o.calculate_center_of_mass(0)
        #
        expected_com = numpy.array([83.286, 14.288, 22.003], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)

    '''

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

