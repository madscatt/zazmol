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

class Test_intg_operate_Move_align(unittest.TestCase): 

    def setUp(self):
        warnings.filterwarnings('ignore')

        self.o1 = system.Molecule(0)
        self.o2 = system.Molecule(1)

        self.o1Sub = system.Molecule(0)
        self.o2Sub = system.Molecule(1)

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

    def test_null_assign(self):
        with self.assertRaises(Exception):
            align_variables = self.o2.align(self.o1, basis_1, basis_2, mode='initialization')

    def test_null_production(self):
        with self.assertRaises(Exception):
            self.o2.align(self.o1, align_variables=align_variables)


    def test_1ATM_against_1ATMRot_subset_full_pdb(self):
        frame = 0
        self.o1.read_pdb(DataPath+'1ATM.pdb')
        self.o2.read_pdb(DataPath+'1ATM.pdb')

        basis_1 = 'resid[i] < 100000'
        basis_2 = 'resid[i] < 100000'

        # Initialization mode
        align_variables = self.o2.align(self.o1, basis_1, basis_2, mode='initialization')

        # Production mode
        self.o2.align(self.o1, align_variables=align_variables)

        expected_coor = self.o1.coor()[frame]
        result_coor = self.o2.coor()[frame]

        self.assert_list_almost_equal(expected_coor, result_coor,places=3)

    def test_2AAD_against_2AADRot_subset_full_pdb(self):
        frame = 0

        self.o1.read_pdb(DataPath+'2AAD.pdb')
        self.o2.read_pdb(DataPath+'2AAD.pdb')

        basis_1 = 'resid[i] < 100000'
        basis_2 = 'resid[i] < 100000'

        # Initialization mode
        align_variables = self.o2.align(self.o1, basis_1, basis_2, mode='initialization')

        # Production mode
        self.o2.align(self.o1, align_variables=align_variables)

        expected_coor = self.o1.coor()[frame]
        result_coor = self.o2.coor()[frame]

        self.assert_list_almost_equal(expected_coor, result_coor,places=3)

    def test_1CRN_against_1CRNRot_subset_full_pdb(self):
        frame = 0
        self.o1.read_pdb(DataPath+'1CRN.pdb')
        self.o2.read_pdb(moduleDataPath+'1CRN-rot.pdb')

        basis_1 = 'resid[i] < 100000'
        basis_2 = 'resid[i] < 100000'

        # Initialization mode
        align_variables = self.o2.align(self.o1, basis_1, basis_2, mode='initialization')

        # Production mode
        self.o2.align(self.o1, align_variables=align_variables)

        expected_coor = self.o1.coor()[frame]
        result_coor = self.o2.coor()[frame]

        self.assert_list_almost_equal(expected_coor, result_coor,places=2)


    def test_1CRN_against_1CRNRotShift_subset_full_pdb(self):
        frame = 0
        self.o1.read_pdb(DataPath+'1CRN.pdb')
        self.o2.read_pdb(moduleDataPath+'1CRN-rot-shift.pdb')

        basis_1 = 'resid[i] < 100000'
        basis_2 = 'resid[i] < 100000'

        # Initialization mode
        align_variables = self.o2.align(self.o1, basis_1, basis_2, mode='initialization')

        # Production mode
        self.o2.align(self.o1, align_variables=align_variables)

        expected_coor = self.o1.coor()[frame]
        result_coor = self.o2.coor()[frame]

        self.assert_list_almost_equal(expected_coor, result_coor,places=2)


    def test_1CRN_pdb(self):
        '''
        read in the same PDB file into two separate molecules
        rotate the second molecule and 
        set a basis then align the second molecule on the first molecule
        '''
        frame = 0

        self.o1.read_pdb(DataPath + '1CRN.pdb')
        self.o2.read_pdb(DataPath + '1CRN.pdb')

        basis_1 = 'name[i] == "CA" and (resid[i] >= 20 and resid[i] <= 31)'
        basis_2 = 'name[i] == "CA" and (resid[i] >= 20 and resid[i] <= 31)'

        # Initialization mode
        align_variables = self.o2.align(self.o1, basis_1, basis_2, mode='initialization')

        axis = 'z'
        theta = numpy.pi / 2.0

        self.o2.rotate(frame, axis, theta)
      
        # Production mode
        self.o2.align(self.o1, align_variables=align_variables)

        expected_com = self.o1.calculate_center_of_mass(frame)
        result_com = self.o2.calculate_center_of_mass(frame)

        self.assert_list_almost_equal(expected_com, result_com, 2)

    def tearDown(self):
        pass

if __name__ == '__main__': 
    main()
