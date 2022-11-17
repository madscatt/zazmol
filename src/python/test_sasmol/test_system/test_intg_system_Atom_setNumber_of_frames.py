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

from unittest import main 
import unittest

import numpy, os, copy

import warnings; warnings.filterwarnings('ignore')

floattype=os.environ['SASMOL_FLOATTYPE']

import sasmol.system as system

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','system')+os.path.sep

class Test_intg_system_Atom_setNumber_of_frames(unittest.TestCase):

   def setUp(self):
      self.o=system.Atom(3,'1CRN-3frames.pdb')

   def test_1CRN_3frames(self):
      '''
	   test a regular pdb file with 3 frame
	   '''
      #
      number_of_frames = 3
      expected_number_of_frames = number_of_frames
      #
      self.o.read_pdb(DataPath+'1CRN-3frames.pdb')
      self.o.setNumber_of_frames(number_of_frames)
      #
      result_number_of_frames = self.o.number_of_frames()
      #
      self.assertEqual(expected_number_of_frames, result_number_of_frames)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   unittest.main() 

