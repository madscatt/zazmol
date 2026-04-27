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

from unittest import main,skipIf 
import unittest
import sasmol.system as system

import numpy, os, copy

import warnings; warnings.filterwarnings('ignore')

floattype=os.environ['SASMOL_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','file_io')+os.path.sep


class Test_intg_file_io_Files_write_pdb(unittest.TestCase):

   def setUp(self):
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


   def test_1ATM(self):
      '''
	   test a 1-atom pdb
	   '''
      #
      frame = 0
      self.o.read_pdb(DataPath+'1ATM.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/1ATM-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1ATM-writepdb-test.pdb')


   def test_1ATM_frame1to2_second_frame(self):
      '''
	   test the second frame a pdb file with 1 atom and 2 frames
	   '''
      #
      frame = 1
      self.o.read_pdb(DataPath+'1ATM-1to2.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/1ATM-1to2-2-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1ATM-1to2-2-writepdb-test.pdb')


   def test_2AAD(self):
      '''
	   test a 2-aa pdb
	   '''
      #
      frame = 0
      self.o.read_pdb(DataPath+'2AAD.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/2AAD-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/2AAD-writepdb-test.pdb')


   def test_2AAD_frames1to3_second_frame(self):
      '''
	   test the second frame a pdb file with 2 amino acids and 3 frames
	   '''
      #
      frame = 2
      self.o.read_pdb(DataPath+'2AAD-1to3.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/2AAD-1to3-2-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/2AAD-1to3-2-writepdb-test.pdb')

   def test_rna_frame1to10_frame_3(self):
      '''
	   test the 10th frame of a pdb file of rna with 10 frames
      '''
      #
      frame = 3
      self.o.read_pdb(DataPath+"rna-1to10.pdb")
      result = self.o.write_pdb(moduleDataPath+'test-results/rna-1to10-3-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/rna-1to10-3-writepdb-test.pdb')


   def test_1CRN(self):
      '''
	   test a small protein (crambin)
      '''
      #
      frame = 0
      self.o.read_pdb(DataPath+"1CRN.pdb")
      result = self.o.write_pdb(moduleDataPath+'test-results/1CRN-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1CRN-writepdb-test.pdb')


   def test_missing_optional_fields_are_not_filled_by_default(self):
      '''
      test write_pdb default behavior does not fill missing optional fields
      '''
      filename = moduleDataPath+'test-results/missing-optional-default.pdb'
      self.o = system.Molecule_Maker(1)
      self.o.setLoc([])
      self.o.setRescode([])
      self.o.setOccupancy([])
      self.o.setBeta([])
      self.o.setSegname([])
      self.o.setElement([])
      self.o.setCharge([])

      result = self.o.write_pdb(filename, 0, 'w')

      with open(filename) as outfile:
         lines = outfile.readlines()
      self.assertEqual(result,1)
      self.assertEqual([line for line in lines if line.startswith('ATOM')], [])
      os.remove(filename)


   def test_missing_optional_fields_can_be_filled_for_write_pdb(self):
      '''
      test write_pdb can fill missing optional fields when explicitly requested
      '''
      filename = moduleDataPath+'test-results/missing-optional-filled.pdb'
      self.o = system.Molecule_Maker(1)
      self.o.setLoc([])
      self.o.setRescode([])
      self.o.setOccupancy([])
      self.o.setBeta([])
      self.o.setSegname([])
      self.o.setElement([])
      self.o.setCharge([])

      result = self.o.write_pdb(filename, 0, 'w', fill_missing_optional=True)

      with open(filename) as outfile:
         lines = outfile.readlines()
      atom_lines = [line for line in lines if line.startswith('ATOM')]
      self.assertEqual(result,1)
      self.assertEqual(len(atom_lines), 1)
      self.assertEqual(atom_lines[0][16], ' ')
      self.assertEqual(atom_lines[0][26], ' ')
      self.assertEqual(atom_lines[0][54:60], '  0.00')
      self.assertEqual(atom_lines[0][60:66], '  0.00')
      self.assertEqual(atom_lines[0][72:76], '    ')
      self.assertEqual(atom_lines[0][76:78], '  ')
      self.assertEqual(atom_lines[0][78:80], '  ')
      os.remove(filename)


   def test_missing_required_fields_raise_when_filling_optional_fields(self):
      '''
      test fill_missing_optional does not fill required PDB fields
      '''
      filename = moduleDataPath+'test-results/missing-required.pdb'
      self.o = system.Molecule_Maker(1)
      self.o.setName([])

      with self.assertRaises(IndexError):
         self.o.write_pdb(filename, 0, 'w', fill_missing_optional=True)
      os.remove(filename)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")   
   def test_1KP8(self):
      '''
	   test a large protein complex (groel)
      '''
      #
      frame = 0
      self.o.read_pdb(DataPath+"1KP8.pdb")
      result = self.o.write_pdb(moduleDataPath+'test-results/1KP8-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1KP8-writepdb-test.pdb')



   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   unittest.main() 

