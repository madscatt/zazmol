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
import sasmol.dcdio as dcdio

import warnings

import os, sys, string, shutil

commonDcdDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','dcd_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','file_io')+os.path.sep


class Test_intg_sasio_Files_close_dcd_write(unittest.TestCase):

   def setUp(self):
      self.o=system.Molecule(0)
      warnings.filterwarnings('ignore')
      
   def test_file_doesnt_exist(self):
      '''
	   test a dcd which doent exist
	   '''
      filename = 'file-notexist.dcd'
      dcdFile = moduleDataPath+'test-results/'+filename
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(dcdFile)
      #sys.stdout = open(stdoutFile,'w')
      with open(stdoutFile,'w') as sys.stdout:
        self.o.close_dcd_write(pf)

      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(dcdFile)



   def test_1ATM(self):
      '''
	   test a dcd with 2 frames based on an 1-atom pdb
	   '''
      #
      filename = '1ATM.dcd'
      dcdFile = commonDcdDataPath+filename
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      with open(stdoutFile,'w') as sys.stdout:
        self.o.close_dcd_write(pf)
        
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)



   def test_2AAD(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb

	   '''
      #
      filename = '2AAD.dcd'
      dcdFile = commonDcdDataPath+filename
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)


   def test_rna_1to10(self):
      '''
	   test a dcd with 10 frames based on a 2-aa pdb

	   '''
      #
      filename = 'rna-1to10.dcd'
      dcdFile = commonDcdDataPath+filename
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)


   @skipIf(os.environ['SASMOL_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_0point8gb(self):
      '''
	   test a dcd of size 0.8gb based on a rna molecule
	   '''
      #
      filename = "rna-0.8g.dcd"
      dcdFile = os.path.join('/tmp/',filename)
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)


   @skipIf(os.environ['SASMOL_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_1point0gb(self):
      '''
	   test a dcd of size 1.0gb based on a rna molecule
	   '''
      #
      filename = "rna-1.0g.dcd"
      dcdFile = os.path.join('/tmp/',filename)
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)



   @skipIf(os.environ['SASMOL_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_2point0gb(self):
      '''
	   test a dcd of size 2.0gb based on a rna molecule
	   '''
      #
      filename = "rna-2.0g.dcd"
      dcdFile = os.path.join('/tmp/',filename)
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)


   
   @skipIf(os.environ['SASMOL_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_3point2gb(self):
      '''
	   test a dcd of size 3.2gb based on a rna molecule
	   '''
      #
      filename = "rna-3.2g.dcd"
      dcdFile = os.path.join('/tmp/',filename)
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)


   
   @skipIf(os.environ['SASMOL_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_6point4gb(self):
      '''
	   test a dcd of size 6.4gb based on a rna molecule
	   '''
      #
      filename = "rna-6.4g.dcd"
      dcdFile = os.path.join('/tmp/',filename)
      tmpDcdFile = moduleDataPath+'test-results/'+filename
      shutil.copy(dcdFile, tmpDcdFile)
      stdoutFile = filename+'.stdiout'
      pf = dcdio.open_dcd_write(tmpDcdFile)
      sys.stdout = open(stdoutFile,'w')
      self.o.close_dcd_write(pf)
      sys.stdout = sys.__stdout__
      code=open(stdoutFile).read().split()[2]
      self.assertEqual(int(code), 0)
      os.remove(stdoutFile)
      os.remove(tmpDcdFile)



   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   unittest.main() 

