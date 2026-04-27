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

import sasmol.system as system


class Test_intg_system_Atom_filename(unittest.TestCase):

    def setUp(self):
        self.o = system.Atom(3, '1CRN-3frames.pdb')

    def test_filename(self):
        expected = '1CRN.pdb'
        self.o.setFilename(expected)
        result = self.o.filename()
        self.assertEqual(expected, result)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
