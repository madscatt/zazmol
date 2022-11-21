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


from unittest import main 
import unittest
import sasmol.file_io as file_io

class Test_intg_sasmol_Files_init(unittest.TestCase):

   def setUp(self):
      pass

   def test_init(self):
      '''
      test initializer of sasmol.Files
      '''
      #
      filename = 'null' #there is really nothing in the constructor
      flag = 0
      o=file_io.Files(filename, flag)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
    unittest.main() 

