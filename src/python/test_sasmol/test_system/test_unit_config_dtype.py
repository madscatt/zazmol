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

import unittest

import numpy

import sasmol.config as config


class Test_config_dtype(unittest.TestCase):

    def test_coordinate_dtype(self):
        '''
        coordinate storage dtype is single precision
        '''
        #
        self.assertEqual(config.COORD_DTYPE, numpy.float32)

    def test_calculation_dtype(self):
        '''
        derived calculation dtype is double precision
        '''
        #
        self.assertEqual(config.CALC_DTYPE, numpy.float64)


if __name__ == '__main__':
    unittest.main()
