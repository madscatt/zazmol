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
from mocker import Mocker, MockerTestCase, ANY, ARGS
import sasmol.system as system

import numpy

import os
floattype=os.environ['SASMOL_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_sascalc_Prop_calccom(MockerTestCase): 

    def setUp(self):
        self.o=system.Molecule(0)

    def calc_exp(self):
        self.o.calculate_center_of_mass(0)
        coor = numpy.array((self.o._coor[0]),floattype)
        return numpy.sqrt(numpy.sum((coor-self.o._com)**2)/len(coor))


    def test_null(self):
        with self.assertRaises(Exception):
            self.o.calculate_radius_of_gyration(0)

    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        self.o.setTotal_mass(0.0)
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calculate_radius_of_gyration(0)
        expected_rg = 0.0
        self.assertAlmostEqual(expected_rg, result_rg)

    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        self.o.setTotal_mass(0.0) 
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calculate_radius_of_gyration(0)
        expected_rg = 2.998744
        self.assertAlmostEqual(expected_rg, result_rg, 5)

    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        self.o.setTotal_mass(0.0) 
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calculate_radius_of_gyration(0)
        expected_rg = self.calc_exp()
        print(result_rg, expected_rg)
        self.assertAlmostEqual(expected_rg, result_rg, 3)

    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        self.o.setTotal_mass(0.0) 
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calculate_radius_of_gyration(0)
        expected_rg = self.calc_exp()
        print(result_rg, expected_rg)
        self.assertAlmostEqual(expected_rg, result_rg, 3)

    @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        self.o.setTotal_mass(0.0) 
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calculate_radius_of_gyration(0)
        expected_rg = self.calc_exp()
        print(result_rg, expected_rg)
        self.assertAlmostEqual(expected_rg, result_rg, 3)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

