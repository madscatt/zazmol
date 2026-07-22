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

from collections import Counter
from unittest import main, skipIf
import unittest
import sasmol.system as system

import numpy, os, copy, tempfile

floattype=os.environ['SASMOL_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','file_io')+os.path.sep

class Test_intg_file_io_Files_read_dcd(unittest.TestCase):
   '''
   Characterize intentionally tolerant PDB reader behavior.

   These tests protect accepted real-world inputs such as blank trailing lines,
   single-frame files without END records, END-separated pseudo-trajectories,
   MODEL/ENDMDL trajectories, non-CHARMM atom names, and non-protein systems.
   '''

   def setUp(self):
      self.o=system.Molecule(0)

   def test_file_doesnt_exist(self):
      '''
      test that missing pdb files raise a Python file error
      '''
      pdbFileName = DataPath+'file-notexist.pdb'

      with self.assertRaises(FileNotFoundError):
         self.o.read_pdb(pdbFileName)

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

   def make_atom_line(self, serial, name, resname, chain, resid, x, segname=None):
      if segname is None:
         segname = chain
      return "%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8.3f%8.3f%8.3f%6s%6s      %-4s%2s%2s\n" % (
         'ATOM', serial, name, ' ', resname, chain, resid, ' ',
         x, x + 1.0, x + 2.0, '1.00', '0.00', segname, name[0], '  ')

   def write_temp_pdb(self, lines):
      handle = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
      try:
         handle.writelines(lines)
      finally:
         handle.close()
      return handle.name

   def test_resname_moltype_classification_does_not_resolve_overlap_as_rna(self):
      lines = [
         self.make_atom_line(1, 'P', 'ALA', 'P', 1, 1.0),
         self.make_atom_line(2, 'P', 'ADE', 'D', 2, 2.0),
         self.make_atom_line(3, 'P', 'GUA', 'D', 3, 3.0),
         self.make_atom_line(4, 'P', 'CYT', 'D', 4, 4.0),
         self.make_atom_line(5, 'P', 'THY', 'D', 5, 5.0),
         self.make_atom_line(6, 'P', 'URA', 'R', 6, 6.0),
         self.make_atom_line(7, 'P', 'A', 'R', 7, 7.0),
         self.make_atom_line(8, 'P', 'DA', 'D', 8, 8.0),
         self.make_atom_line(9, 'O', 'TIP3', 'W', 9, 9.0),
         self.make_atom_line(10, 'C1', 'LIG', 'L', 10, 10.0),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(
         list(self.o.moltype()),
         ['protein', 'dna', 'dna', 'dna', 'dna',
          'nucleic', 'nucleic', 'dna', 'water', 'other'])
      self.assertEqual(Counter(self.o.moltype()),
                       Counter({'dna': 5, 'nucleic': 2,
                                'protein': 1, 'water': 1, 'other': 1}))

   def test_o2prime_evidence_promotes_overlap_segment_to_rna(self):
      lines = [
         self.make_atom_line(1, "O2'", 'ADE', 'R', 1, 1.0),
         self.make_atom_line(2, 'P', 'GUA', 'R', 2, 2.0),
         self.make_atom_line(3, 'P', 'CYT', 'R', 3, 3.0),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(list(self.o.moltype()), ['rna', 'rna', 'rna'])

   def test_o2star_alias_promotes_overlap_segment_to_rna(self):
      lines = [
         self.make_atom_line(1, 'O2*', 'ADE', 'R', 1, 1.0),
         self.make_atom_line(2, 'P', 'GUA', 'R', 2, 2.0),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(list(self.o.moltype()), ['rna', 'rna'])

   def test_ambiguous_thymine_does_not_conflict_with_o2prime_rna(self):
      lines = [
         self.make_atom_line(1, "O2'", 'ADE', 'X', 1, 1.0),
         self.make_atom_line(2, 'P', 'THY', 'X', 2, 2.0),
         self.make_atom_line(3, 'P', 'GUA', 'X', 3, 3.0),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(list(self.o.moltype()), ['rna', 'rna', 'rna'])
      report = self.o.moltype_by_segname_report()
      self.assertEqual(report['overall_status'], 'clean')
      self.assertEqual(report['segments']['X']['identity'], 'resolved_rna')

   def test_pdbscan_blank_segnames_use_chain_for_o2prime_inference(self):
      lines = [
         self.make_atom_line(1, "O2'", 'ADE', 'R', 1, 1.0, segname=''),
         self.make_atom_line(2, 'P', 'GUA', 'R', 2, 2.0, segname=''),
         self.make_atom_line(3, 'P', 'ADE', 'D', 3, 3.0, segname=''),
         self.make_atom_line(4, 'P', 'GUA', 'D', 4, 4.0, segname=''),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file, pdbscan=True)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(list(self.o.segname()), ['', '', '', ''])
      self.assertEqual(list(self.o.moltype()), ['rna', 'rna', 'nucleic', 'nucleic'])
      report = self.o.moltype_by_segname_report()
      self.assertEqual(report['segments']['chain:R']['grouping_source'], 'chain')
      self.assertEqual(report['segments']['chain:R']['chain'], 'R')
      self.assertEqual(report['segments']['chain:D']['status'], 'ambiguous')

   def test_pdbscan_chain_fallback_isolates_cross_chain_o2prime_conflict(self):
      lines = [
         self.make_atom_line(1, "O2'", 'ADE', 'R', 1, 1.0, segname=''),
         self.make_atom_line(2, 'P', 'GUA', 'R', 2, 2.0, segname=''),
         self.make_atom_line(3, 'P', 'THY', 'D', 3, 3.0, segname=''),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file, pdbscan=True)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(list(self.o.moltype()), ['rna', 'rna', 'nucleic'])
      report = self.o.moltype_by_segname_report()
      self.assertEqual(report['overall_status'], 'ambiguous')
      self.assertEqual(report['segments']['chain:R']['status'], 'resolved_rna')
      self.assertEqual(report['segments']['chain:D']['status'], 'ambiguous')

   def test_pdbscan_blank_segname_and_blank_chain_does_not_propagate_o2prime(self):
      lines = [
         self.make_atom_line(1, "O2'", 'ADE', '', 1, 1.0, segname=''),
         self.make_atom_line(2, 'P', 'GUA', '', 2, 2.0, segname=''),
         'END\n']
      pdb_file = self.write_temp_pdb(lines)
      try:
         self.o.read_pdb(pdb_file, pdbscan=True)
      finally:
         os.unlink(pdb_file)

      self.assertEqual(list(self.o.moltype()), ['rna', 'nucleic'])
      report = self.o.moltype_by_segname_report()
      self.assertEqual(report['segments']['unassigned:1']['grouping_source'], 'none')

   def test_1ATM_one_frame(self):
      '''
	   test a pdb file with 1 atom and 1 frame
	   '''
      #
      self.o.read_pdb(DataPath+'1ATM.pdb')
      result_coor = self.o.coor()
      #print('\nresult_coor \n',result_coor)
      #
      expected_coor = numpy.array([[[73.944, 41.799, 41.652]]],floattype)
      #print('\nexpected_coor \n',expected_coor)
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)

   def test_vmd_hex_ids_require_explicit_read_mode(self):
      '''
      Hex-like VMD overflow IDs are not accepted by default.
      '''
      pdb_file = self.write_temp_pdb([
         self.make_atom_line('186A0', 'CA', 'GLY', 'A', '271A', 1.0),
         'END\n'])

      try:
         with self.assertRaises(Exception):
            self.o.read_pdb(pdb_file)
      finally:
         os.unlink(pdb_file)

   def test_vmd_hex_read_mode_decodes_unambiguous_atom_and_resid(self):
      '''
      Explicit VMD mode decodes A-F atom serial and resid overflow fields.
      '''
      pdb_file = self.write_temp_pdb([
         self.make_atom_line('186A0', 'CA', 'GLY', 'A', '271A', 1.0),
         self.make_atom_line('186A1', 'CB', 'GLY', 'A', '271A', 2.0),
         'END\n'])

      try:
         self.o.read_pdb(pdb_file, pdb_id_mode='vmd_hex')
         self.assert_list_almost_equal(self.o.index(), [1, 2])
         self.assert_list_almost_equal(self.o.original_index(), [100000, 100001])
         self.assert_list_almost_equal(self.o.resid(), [10010, 10010])
         self.assert_list_almost_equal(self.o.original_resid(), [10010, 10010])
      finally:
         os.unlink(pdb_file)

   def test_vmd_hex_read_mode_decodes_numeric_resid_overflow_with_context(self):
      '''
      VMD resid 10000 is raw 2710, so context is needed to avoid misreading it.
      '''
      pdb_file = self.write_temp_pdb([
         self.make_atom_line('99999', 'CA', 'GLY', 'A', '9999', 1.0),
         self.make_atom_line('186A0', 'CB', 'GLY', 'A', '2710', 2.0),
         'END\n'])

      try:
         self.o.read_pdb(pdb_file, pdb_id_mode='vmd_hex')
         self.assert_list_almost_equal(self.o.original_index(), [99999, 100000])
         self.assert_list_almost_equal(self.o.resid(), [9999, 10000])
      finally:
         os.unlink(pdb_file)

   def test_vmd_hex_read_mode_rejects_ambiguous_numeric_starting_resid(self):
      '''
      A file that starts at raw resid 2710 may mean decimal 2710 or hex 10000.
      '''
      pdb_file = self.write_temp_pdb([
         self.make_atom_line('186A0', 'CA', 'GLY', 'A', '2710', 1.0),
         'END\n'])

      try:
         with self.assertRaises(Exception):
            self.o.read_pdb(pdb_file, pdb_id_mode='vmd_hex')
      finally:
         os.unlink(pdb_file)

   def test_vmd_hex_read_mode_can_force_ambiguous_numeric_fields_to_hex(self):
      '''
      Callers that know a file is VMD-hex encoded can opt into hex ambiguity.
      '''
      pdb_file = self.write_temp_pdb([
         self.make_atom_line('20000', 'CA', 'GLY', 'A', '2710', 1.0),
         'END\n'])

      try:
         self.o.read_pdb(
            pdb_file, pdb_id_mode='vmd_hex', ambiguous_vmd_hex='hex')
         self.assert_list_almost_equal(self.o.original_index(), [131072])
         self.assert_list_almost_equal(self.o.resid(), [10000])
      finally:
         os.unlink(pdb_file)

   def test_vmd_hex_read_mode_decodes_pdbscan_conect(self):
      '''
      pdbscan CONECT parsing follows the same explicit VMD ID mode.
      '''
      pdb_file = self.write_temp_pdb([
         self.make_atom_line('186A0', 'CA', 'GLY', 'A', '271A', 1.0),
         self.make_atom_line('186A1', 'CB', 'GLY', 'A', '271A', 2.0),
         'CONECT186A0186A1\n',
         'END\n'])

      try:
         self.o.read_pdb(pdb_file, pdbscan=True, pdb_id_mode='vmd_hex')
         self.assertEqual(self.o.conect(), {100000: [100001]})
      finally:
         os.unlink(pdb_file)


   def test_1ATM_two_frames(self):
      '''
	   test a pdb file with 1 atom and 2 frames
	   '''
      #
      self.o.read_pdb(DataPath+'1ATM-1to2.pdb')
      result_coor = self.o.coor()
      #print('\nresult_coor \n',result_coor)
      #
      expected_coor = numpy.array([[[76.944, 41.799, 41.652]],[[73.944, 38.799, 41.652]]],floattype)
      #print('\nexpected_coor \n',expected_coor)
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)


   def test_2AAD_one_frame(self):
      '''
	   test a pdb file with 2 amino acids and 1 frame
	   '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      result_coor = self.o.coor()
      #print('\nresult_coor \n',result_coor)
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]]],floattype)
      #print('\nexpected_coor \n',expected_coor)
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)


   def test_2AAD_three_frames_separatedby_END(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by END)
	   '''
      #
      self.o.read_pdb(moduleDataPath+'2AAD-1to3-END.pdb')
      result_coor = self.o.coor()
      #print('\nresult_coor \n',result_coor)
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]],\
      [[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]],\
      [[  73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)
      #print('\nexpected_coor \n',expected_coor)
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)
 

   def test_2AAD_three_frames_separatedby_END_all_properties(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by END) for all properties
	   '''
      #
      self.o.read_pdb(moduleDataPath+'2AAD-1to3-END.pdb')
      natoms = self.o.natoms()
      self.assertEqual(natoms,15)
      self.assertEqual(self.o.moltype(),['protein']*natoms)
      self.assertEqual(self.o.number_of_frames(),3)
      self.assertEqual(self.o.atom(),['ATOM']*natoms)
      self.assertEqual(self.o.name(),['N','CA','C','O','CB','CG1','CG2','CD1','N','CA','C','O','CB','OG1','CG2'])
      self.assert_list_almost_equal(self.o.index(),list(range(1,natoms+1)))
      self.assertEqual(self.o.loc(),[' ']*natoms)
      self.assertEqual(self.o.resname(),['ILE']*8+['THR']*7)
      self.assertEqual(self.o.chain(),['N']*natoms)
      self.assert_list_almost_equal(self.o.resid(),[515]*8+[516]*7)
      self.assertEqual(self.o.rescode(),[' ']*natoms)
      self.assertEqual(self.o.occupancy(),['1.00']*natoms)
      self.assertEqual(self.o.beta(),['36.37', '36.23', '36.32', '36.04', '36.69', '38.12', '34.42', '39.85', '35.01', '35.51', '38.09', '36.94', '36.74', '37.19', '34.44'])
      self.assertEqual(self.o.segname(),['N']*natoms)
      self.assertEqual(self.o.element(),['N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'O', 'C'])
      self.assertEqual(self.o.charge(),['  ']*natoms)


   def test_2AAD_three_frames_separatedby_END_wrong_number_atoms(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-END_wrong_number_atoms.pdb')



   def test_2AAD_three_frames_separatedby_MODEL(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL.pdb')
      result_coor = self.o.coor()
      #print('\nresult_coor \n',result_coor)
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]],\
      [[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]],\
      [[  73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)
      #print('\nexpected_coor \n',expected_coor)
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)
 

   def test_2AAD_three_frames_separatedby_MODEL_wrong_number_atoms(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):      
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL_wrong_number_atoms.pdb')



   def test_2AAD_three_frames_separatedby_MODEL_wrongnumber_mix_END(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL_wrongnumber_mix_END.pdb')



   def test_2AAD_three_frames_separatedby_MODEL_mix_END_noterminating(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):      
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL_mix_END_noterminating.pdb')


   
   def test_rna_frame1to10_frame_3(self):
      '''
	   test a pdb file of rna with 10 frames
      '''
      #
      self.o.read_pdb(DataPath+"rna-1to10.pdb")
      result_coor = self.o.coor()
      #print('\nresult_coor \n',result_coor[2][10627])
      #
      self.assertEqual(len(result_coor),10)
      self.assertEqual(len(result_coor[3]),10632)
      expected_coor_sample = numpy.array([-5.564, 20.324, 26.185],floattype) #atom 10627 of frame 3
      self.assert_list_almost_equal(result_coor[2][10627],expected_coor_sample,3)


   def test_1PSI(self):
      '''
	   test a pdb file without ENDMDL
	   '''
      #
      with self.assertRaises(Exception):
         self.o.read_pdb(DataPath+'1PSI.pdb')

   def test_blanklines(self):
      '''
	   test a pdb file ending with blank lines
	   '''
      #
      self.o.read_pdb(DataPath+'dimcd_fixed_atoms.pdb')
      expected_coor_sample = numpy.array([65.124,  35.624,  50.733],floattype)
      result_coor = self.o.coor()
      self.assert_list_almost_equal(result_coor[0][10],expected_coor_sample,3)


   def test_1AA_NoEND(self):
      '''
	   test a 1AA pdb file with 1frame and without END statement
	   '''
      #
      self.o.read_pdb(moduleDataPath+'1AA-NoEND.pdb')
      result_coor = self.o.coor()
      result_sum_coor = sum(sum(sum((result_coor))))
      #
      expected_coor = numpy.array([[-21.525, -67.562,  86.759], [-22.003, -68.460, 86.892],[-21.905, -66.929,  87.525],[-20.492, -67.726, 86.876],[-21.725, -66.910, 85.457],[-21.476, -67.600, 84.661],[-21.157, -65.997, 85.450],[-23.103, -66.411, 85.215],[-23.249, -65.504, 84.385]],floattype)
      expected_sum_coor = sum(sum(expected_coor))
      #
      self.assertAlmostEqual(result_sum_coor,expected_sum_coor,3)

   def test_cleaned_up_package_rna(self):
      '''
	   test a pdb file of rna with 1250 frames of size 1.0g
      '''
      #
      self.o.read_pdb(moduleDataPath+"new_package_rna.pdb")
      result_coor = self.o.coor()
      #print('\nlength of result_coor \n',len(result_coor[0]))
      #print('\nresult_coor \n',result_coor)
      #
      self.assertEqual(len(result_coor[0]),3719)
      expected_coor_sample = numpy.array([-12.872, 13.360, -153.873],floattype) #atom 10627 of frame 3
      self.assert_list_almost_equal(result_coor[0][299],expected_coor_sample,3)

   def test_problem_pdb(self):
      '''
	   test a pdb file with non-charmm atom names
      '''
      #
      #print('ZHL')
      self.o.read_pdb(moduleDataPath+"nef_nohis.pdb")
      #print(self.o.name())

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   unittest.main() 
