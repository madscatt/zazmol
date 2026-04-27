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
Integration test for sasio.subset.Mask.merge_two_molecules

contract:

merge a non-existing pdb with another non-existing pdb
merge a non-existing pdb with a pdb of 1 atom

merge a pdb of 1 atom with a non-existing pdb
merge a pdb of 1 atom with a problem pdb (1PSI)
merge a pdb of 1 atom with itself
merge a pdb of 1 atom with 1 atom and 2 frames
merge a pdb of 1 atom with a pdb of 2 residues

merge a pdb of 2 residues with a non-existing pdb
merge a pdb of 2 residues with a problem pdb (1PSI)
merge a pdb of 2 residues with itself
merge a pdb of 2 residues with a pdb of 2 residues and 3 frames
merge a pdb of 2 residues with small protein


merge a small protein (crambin) with a non-existing pdb
merge a small protein (crambin) with a problem pdb (1PSI) 
merge a small protein (crambin) with itself
merge a small protein (crambin) with a large protein (groel) (Skipped as SASMOL_LARGETEST)
merge a small protein (crambin) with an rna

merge a large protein (groel) with a non-existing pdb  (Skipped as SASMOL_LARGETEST)
merge a large protein (groel) with a problem pdb (1PSI)  (Skipped as SASMOL_LARGETEST)
merge a large protein (groel) with itself (Skipped as SASMOL_LARGETEST)
merge a large protein (groel) with rna (Skipped as SASMOL_LARGETEST)

merge an rna with a non-existing pdb
merge an rna with a problem pdb (1PSI) 
merge an rna with itself
merge an rna with a large protein complex (groel) (Skipped as SASMOL_LARGETEST)

merge a problem pdb (1PSI) with a non-existing pdb
merge a problem pdb (1PSI) with itself
merge a problem pdb (1PSI) with a large protein complex (groel) (Skipped as SASMOL_LARGETEST)
"""

from unittest import main,skipIf 
import unittest

import sasmol.system as system
import sasmol.config as config
import sasmol.subset as subset

import numpy

import os

PdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
#print('PdbDataPath = ', PdbDataPath)
#print('PdbDataPath = ', PdbDataPath)
#print('PdbDataPath = ', PdbDataPath)

class Test_subset_Mask_merge_two_molecules(unittest.TestCase): 
 

   def setUp(self):
      self.o1=system.Molecule(0)
      self.o2=system.Molecule(1)
      self.o3=system.Molecule(2)

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


   def assert_pdb(self, o1, o2, o3):
      frame = 0
      try:
         self.assertEqual(o1.atom()+o2.atom(), o3.atom())
         self.assert_list_almost_equal(list(o3.index()), list(range(1,len(list(o3.index()))+1)))
         self.assertEqual(o1.name()+o2.name(), o3.name())
         self.assertEqual(o1.loc()+o2.loc(), o3.loc())
         self.assertEqual(o1.resname()+o2.resname(), o3.resname())
         self.assertEqual(o1.chain()+o2.chain(), o3.chain())
         self.assertEqual(list(o1.resid())+list(o2.resid()), list(o3.resid()))
         self.assertEqual(o1.rescode()+o2.rescode(), o3.rescode())
         self.assert_list_almost_equal(list(o1.coor()[frame,:,0])+list(o2.coor()[frame,:,0]), list(o3.coor()[frame,:,0]))  
         self.assert_list_almost_equal(list(o1.coor()[frame,:,1])+list(o2.coor()[frame,:,1]), list(o3.coor()[frame,:,1]))  
         self.assert_list_almost_equal(list(o1.coor()[frame,:,2])+list(o2.coor()[frame,:,2]), list(o3.coor()[frame,:,2]))  
         self.assertEqual(o1.occupancy()+o2.occupancy(), o3.occupancy())
         self.assertEqual(o1.beta()+o2.beta(), o3.beta())
         self.assertEqual(o1.segname()+o2.segname(), o3.segname())
         self.assertEqual(o1.element()+o2.element(), o3.element())
         self.assertEqual(o1.charge()+o2.charge(), o3.charge())
         self.assertEqual(o1.moltype()+o2.moltype(), o3.moltype())
      except:
         raise Exception

   def test_null_1(self):
      '''
      merge a non-existing pdb with another non-existing pdb
      '''
      #
      try:
         self.o1.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      try:
         self.o2.read_pdb(PdbDataPath+'YYY.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('mol1 has no atoms' in error[0])


   def test_null_2(self):
      '''
      merge a non-existing pdb with a pdb of 1 atom
      '''
      #
      try:
         self.o1.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      self.o2.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('mol1 has no atoms' in error[0])


   def test_1ATM_1(self):
      '''
      merge a pdb of 1 atom with a non-existing pdb
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1ATM.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1ATM_2(self):
      '''
      merge a pdb of 1 atom with a problem pdb (1PSI)
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1ATM.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1ATM_3(self):
      '''
      merge a pdb of 1 atom with itself
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1ATM.pdb')
      self.o2.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1ATM_4(self):
      '''
      merge a pdb of 1 atom with 1 atom and 2 frames
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1ATM.pdb')
      self.o2.read_pdb(PdbDataPath+'1ATM-1to2.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1ATM_5(self):
      '''
      merge a pdb of 1 atom with a pdb of 2 residues
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1ATM.pdb')
      self.o2.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_2AAD_1(self):
      '''
      merge a pdb of 2 residues with a non-existing pdb
      '''
      #
      self.o1.read_pdb(PdbDataPath+'2AAD.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_2AAD_2(self):
      '''
      merge a pdb of 2 residues with a problem pdb (1PSI)
      '''
      #
      self.o1.read_pdb(PdbDataPath+'2AAD.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_2AAD_3(self):
      '''
      merge a pdb of 2 residues with itself
      '''
      #
      self.o1.read_pdb(PdbDataPath+'2AAD.pdb')
      self.o2.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_2AAD_4(self):
      '''
      merge a pdb of 2 residues with a pdb of 2 residues and 3 frames
      '''
      #
      self.o1.read_pdb(PdbDataPath+'2AAD.pdb')
      self.o2.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_2AAD_5(self):
      '''
      merge a pdb of 2 residues with small protein
      '''
      #
      self.o1.read_pdb(PdbDataPath+'2AAD.pdb')
      self.o2.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1CRN_1(self):
      '''
      merge a small protein (crambin) with a non-existing pdb
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1CRN.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1CRN_2(self):
      '''
      merge a small protein (crambin) with a problem pdb (1PSI) 
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1CRN.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1CRN_3(self):  
      '''
      merge a small protein (crambin) with itself
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1CRN.pdb')
      self.o2.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1CRN_4(self):
      '''
      merge a small protein (crambin) with a large protein (groel)
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1CRN.pdb')
      self.o2.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1CRN_5(self):
      '''
      merge a small protein (crambin) with an rna
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1CRN.pdb')
      self.o2.read_pdb(PdbDataPath+'rna.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_1(self):
      '''
      merge a large protein (groel) with a non-existing pdb
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1KP8.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_2(self):
      '''
      merge a large protein (groel) with a problem pdb (1PSI) 
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1KP8.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_3(self):
      '''
      merge a large protein (groel) with itself
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1KP8.pdb')
      self.o2.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_4(self):
      '''
      merge a large protein (groel) with rna
      '''
      #
      self.o1.read_pdb(PdbDataPath+'1KP8.pdb')
      self.o2.read_pdb(PdbDataPath+'rna.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_rna_1(self):
      '''
      merge an rna with a non-existing pdb
      '''
      #
      self.o1.read_pdb(PdbDataPath+'rna.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_rna_2(self):
      '''
      merge an rna with a problem pdb (1PSI) 
      '''
      #
      self.o1.read_pdb(PdbDataPath+'rna.pdb')
      try:
         self.o2.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      with self.assertRaises(Exception):
         self.assert_pdb(self.o1, self.o2, self.o3)


   def test_rna_3(self):
      '''
      merge an rna with itself
      '''
      #
      self.o1.read_pdb(PdbDataPath+'rna.pdb')
      self.o2.read_pdb(PdbDataPath+'rna.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_rna_4(self):
      '''
      merge an rna with a large protein complex (groel)
      '''
      #
      self.o1.read_pdb(PdbDataPath+'rna.pdb')
      self.o2.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      expecting_error = False
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o1, self.o2, self.o3)


   def test_1PSI_1(self):
      '''
      merge a problem pdb (1PSI) with a non-existing pdb
      '''
      #
      try:
         self.o1.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      try:
         self.o2.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('mol1 has no atoms' in error[0])



   def test_1PSI_2(self):
      '''
      merge a problem pdb (1PSI) with itself
      '''
      #
      try:
         self.o1.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      try:
         self.o2.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('mol1 has no atoms' in error[0])


   @skipIf(os.environ['SASMOL_LARGETEST']=='n',"I am not testing large files")
   def test_1PSI_3(self):
      '''
      merge a problem pdb (1PSI) with a large protein complex (groel)
      '''
      #
      try:
         self.o1.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      self.o2.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('mol1 has no atoms' in error[0])


   def test_molecule_maker_objects(self):
      '''
      merge two Molecule_Maker objects
      '''
      #
      self.o1 = system.Molecule_Maker(2, name=['N', 'CA'], resid=[10, 10])
      self.o2 = system.Molecule_Maker(2, name=['C', 'O'], resid=[20, 20])
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(error, [])
      self.assertEqual(self.o3.natoms(), 4)
      self.assertEqual(self.o3.name(), ['N', 'CA', 'C', 'O'])
      self.assertEqual(list(self.o3.index()), [1, 2, 3, 4])
      self.assertEqual(list(self.o3.resid()), [10, 10, 20, 20])
      self.assertEqual(list(self.o3.original_index()), [1, 2, 1, 2])
      self.assertEqual(list(self.o3.original_resid()), [10, 10, 20, 20])
      self.assertEqual(self.o3.coor().shape, (1, 4, 3))


   def test_molecule_maker_coordinate_dtype(self):
      '''
      merge stores coordinates as single precision
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(error, [])
      self.assertEqual(self.o1.coor().dtype, config.COORD_DTYPE)
      self.assertEqual(self.o2.coor().dtype, config.COORD_DTYPE)
      self.assertEqual(self.o3.coor().dtype, config.COORD_DTYPE)


   def test_missing_descriptor_is_skipped(self):
      '''
      merge skips descriptors that are not present in both molecules
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      del self.o2._element
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(error, [])
      self.assertFalse(hasattr(self.o3, '_element'))
      self.assertEqual(self.o3.natoms(), 2)
      self.assertEqual(list(self.o3.index()), [1, 2])


   def test_missing_descriptor_report_is_optional(self):
      '''
      report skipped descriptors when requested by the caller
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      del self.o2._element
      #
      error = self.o3.merge_two_molecules(
          self.o1, self.o2, report_missing_descriptors=True)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('skipped descriptor _element' in error[0])
      self.assertTrue('mol2' in error[0])
      self.assertFalse(hasattr(self.o3, '_element'))
      self.assertEqual(self.o3.natoms(), 2)


   def test_missing_core_descriptor_is_skipped(self):
      '''
      merge keeps available descriptors when another descriptor is absent
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      del self.o1._name
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(error, [])
      self.assertFalse(hasattr(self.o3, '_name'))
      self.assertEqual(self.o3.atom(), ['ATOM', 'ATOM'])
      self.assertEqual(self.o3.natoms(), 2)


   def test_missing_coordinates_current_behavior(self):
      '''
      characterize current behavior when coordinates are missing
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      self.o2.setCoor(None)
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('coordinates' in error[0])


   def test_coordinate_atom_count_mismatch_fails(self):
      '''
      merge fails when coordinate atom count disagrees with natoms
      '''
      #
      self.o1 = system.Molecule_Maker(2)
      self.o2 = system.Molecule_Maker(1)
      self.o1.setCoor(numpy.zeros((1, 1, 3), config.COORD_DTYPE))
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(len(error), 1)
      self.assertTrue('mol1 coordinate atom count' in error[0])


   def test_optional_descriptors_merge_when_present(self):
      '''
      merge optional force-field descriptors when both molecules define them
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      self.o1._charmm_type = ['C1']
      self.o2._charmm_type = ['C2']
      self.o1.setAtom_charge([0.1])
      self.o2.setAtom_charge([0.2])
      self.o1.setAtom_vdw([1.5])
      self.o2.setAtom_vdw([1.6])
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(error, [])
      self.assertEqual(list(self.o3._charmm_type), ['C1', 'C2'])
      self.assert_list_almost_equal(list(self.o3.atom_charge()), [0.1, 0.2])
      self.assert_list_almost_equal(list(self.o3.atom_vdw()), [1.5, 1.6])


   def test_connectivity_is_copied_without_aliasing_mol1(self):
      '''
      merge copies connectivity from mol1 without sharing its dictionary
      '''
      #
      self.o1 = system.Molecule_Maker(1)
      self.o2 = system.Molecule_Maker(1)
      self.o1._conect = {1: [2]}
      self.o2._conect = {2: [1]}
      #
      error = self.o3.merge_two_molecules(self.o1, self.o2)
      #
      self.assertEqual(error, [])
      self.assertEqual(self.o3.conect(), {1: [2], 2: [1]})
      self.assertIsNot(self.o3.conect(), self.o1._conect)



   def tearDown(self):
      pass


if __name__ == '__main__': 
   main() 
