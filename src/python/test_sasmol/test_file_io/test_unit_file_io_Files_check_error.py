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
from unittest import mock

from sasmol.pdb_io import PDBElementResolutionError
import sasmol.system as system

import os


class Test_unit_sasio_Files_print_error(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def test_no_error(self):
        '''
        test for no error
        '''
        error = []
        self.o.check_error(error)

    def test_error(self):
        '''
        test for error
        '''
        error = ['wrong']
        with self.assertRaises(PDBElementResolutionError) as captured_error:
            self.o.check_error(error)
        self.assertEqual(str(captured_error.exception), 'wrong')

    def test_read_pdb_logs_and_reraises_element_resolution_error(self):
        '''
        test for caller-boundary logging during PDB element resolution failure
        '''
        data_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..',
            'data',
            'pdb_common'
        ) + os.path.sep

        with mock.patch.object(
            self.o,
            'element_filter',
            side_effect=PDBElementResolutionError('wrong')
        ):
            with self.assertLogs('sasmol.pdb_io', level='ERROR') as captured_logs:
                with self.assertRaises(PDBElementResolutionError):
                    self.o.read_pdb(data_path + '1ATM.pdb')

        self.assertTrue(any('Unable to resolve PDB element names while reading' in line
                            for line in captured_logs.output))

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
