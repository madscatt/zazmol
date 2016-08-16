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

"""
Integration test for sasio.subset.Mask.get_coor_using_mask

contract:

null test by providing the empty mask

negative test by providing the wrong mask

1-atom/1-frame, mask no atom
1-atom/1-frame, mask 1 atom
1-atom/2-frames, mask frame-0(1 atom)
1-atom/2-frames, mask frame-1(1 atom)

2-aa/1-frames, mask frame-0(no atom)
2-aa/1-frames, mask frame-0(8 atoms)
2-aa/1-frames, mask frame-0(8 atoms)
2-aa/3-frames, mask frame-0(no atom)
2-aa/3-frames, mask frame-0(8 atoms)
2-aa/3-frames, mask frame-0(all atoms)
2-aa/3-frames, mask frame-1(8 atoms)
2-aa/3-frames, mask frame-2(8 atoms)

rna/1-frame, mask frame-0(no atom)
rna/1-frame, mask frame-0(248 atoms)
rna/1-frame, mask frame-0(all atoms)

small protein (crambin)/1-frame, mask frame-0(no atom)
small protein (crambin)/1-frame, mask frame-0(46 atoms)
small protein (crambin)/1-frame, mask frame-0(all atoms)

large protein (groel)/1-frame, mask frame-0(no atom) (Skipped as SASMOL_LARGETEST)
large protein (groel)/1-frame, mask frame-0(7350 atoms) (Skipped as SASMOL_LARGETEST)
large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASMOL_LARGETEST)
large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASMOL_LARGETEST)

a bad pdb, assertRaises
"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.system as system
import sasmol.subset as subset 
import numpy

import os

PdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_subset_Mask_get_coor_using_mask(MockerTestCase): 
 

   def setUp(self):
      self.o=system.Molecule(0)
      self.o_result=system.Molecule(1)
      self.o_expected=system.Molecule(2)


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



   def test_null(self):
      '''
      null
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      frame = 0
      #
      expecting_error = False
      mask = [0]
      #
      with self.assertRaises(Exception):
         error, result_coor = self.o.get_coor_using_mask(frame, mask)


   def test_wrong(self):
      '''
      negative test by providing the wrong mask
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      frame = 0
      #
      expecting_error = False
      mask = [1,2]
      #
      with self.assertRaises(Exception):
         error, result_coor = self.o.get_coor_using_mask(frame, mask)


   def test_1ATM_mask_none(self):
      '''
	   1-atom/1-frame, mask no atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      frame = 0
      #
      expecting_error = False 
      mask = [0]
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_1ATM_mask_one(self):
      '''
	   1-atom/1-frame, mask 1 atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      frame = 0
      #
      expecting_error = False 
      mask = [1]
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_1ATM_2frames_mask_frame_0(self):
      '''
	   1-atom/2-frames, mask frame-0(1 atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM-1to2.pdb')
      #
      frame = 0
      #
      expecting_error = False 
      mask = [0]
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_1ATM_2frames_mask_frame1(self):
      '''
	   1-atom/2-frames, mask frame-1(1 atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM-1to2.pdb')
      #
      frame = 1
      #
      expecting_error = False 
      mask = [1]
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_1frame_mask_frame0_atomNone(self):
      '''
	   2-aa/1-frames, mask frame-0(no atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B" or resid[i]==12'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_1frame_mask_frame0_8atoms(self):
      '''
	   2-aa/1-frames, mask frame-0(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_1frame_mask_frame0_atoms_all(self):
      '''
	   2-aa/1-frames, mask frame-0(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_3frames_mask_frame0_atoms_none(self):
      '''
	   2-aa/3-frames, mask frame-0(no atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_3frames_mask_frame0_8atoms(self):
      '''
	   2-aa/3-frames, mask frame-0(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_3frames_mask_frame0_allAtoms(self):
      '''
	   2-aa/3-frames, mask frame-0(all atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]!=5'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_3frames_mask_frame1_8atoms(self):
      '''
	   2-aa/3-frames, mask frame-1(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 1
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_2AAD_3frames_mask_frame2_8atoms(self):
      '''
	   2-aa/3-frames, mask frame-2(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 2
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_rna_1(self):
      '''
      rna/1-frame, mask frame-0(no atom)
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_rna_2(self):
      '''
      rna/1-frame, mask frame-0(248 atoms)
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]==20'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_rna_3(self):
      '''
      rna/1-frame, mask frame-0(all atoms)
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]>0'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_1CRN_1(self):
      '''
      small protein (crambin)/1-frame, mask frame-0(no atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)
      

   def test_1CRN_2(self):
      '''
      small protein (crambin)/1-frame, mask frame-0(46 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="N"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)



   def test_1CRN_3(self):
      '''
      small protein (crambin)/1-frame, mask frame-0(all atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)



   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_1(self):
      '''
      large protein (groel)/1-frame, mask frame-0(no atom) (Skipped as SASMOL_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_2(self):
      '''
      large protein (groel)/1-frame, mask frame-0(7350 atoms) (Skipped as SASMOL_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="N"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_3(self):
      '''
      large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASMOL_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      expected_coor = [[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]]
      #
      error, result_coor = self.o.get_coor_using_mask(frame, mask)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_coor, result_coor)


   def test_1PSI(self):
      '''
      a bad pdb, assertRaises
      '''
      try:
         self.o.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      frame = 0
      #
      expecting_error = False 
      mask = [0]*3718   
      #
      with self.assertRaises(Exception):
         error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)



   def tearDown(self):
      pass


if __name__ == '__main__': 
   main() 
