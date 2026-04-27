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

import os
import unittest

import numpy

import sasmol.config as config
import sasmol.system as system


PDB_DATA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '..', 'data', 'pdb_common') + os.path.sep
DCD_DATA_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '..', 'data', 'dcd_common') + os.path.sep


class Test_coordinate_dtype_characterization(unittest.TestCase):

    def test_molecule_maker_coordinates_are_policy_dtype(self):
        '''
        Molecule_Maker currently stores coordinates as single precision
        '''
        #
        molecule = system.Molecule_Maker(1)
        #
        self.assertEqual(molecule.coor().dtype, config.COORD_DTYPE)

    def test_read_pdb_coordinates_are_policy_dtype(self):
        '''
        read_pdb stores coordinates as single precision
        '''
        #
        molecule = system.Molecule(0)
        molecule.read_pdb(PDB_DATA_PATH + '1ATM.pdb')
        #
        self.assertEqual(molecule.coor().dtype, config.COORD_DTYPE)

    def test_read_dcd_coordinates_are_policy_dtype(self):
        '''
        read_dcd stores coordinates as single precision
        '''
        #
        molecule = system.Molecule(0)
        molecule.read_dcd(DCD_DATA_PATH + '1ATM.dcd')
        #
        self.assertEqual(molecule.coor().dtype, config.COORD_DTYPE)

    def test_read_single_dcd_step_coordinates_are_policy_dtype(self):
        '''
        read_single_dcd_step currently stores coordinates as single precision
        '''
        #
        molecule = system.Molecule(0)
        molecule.read_single_dcd_step(DCD_DATA_PATH + '1ATM.dcd', 1)
        #
        self.assertEqual(molecule.coor().dtype, config.COORD_DTYPE)

    def test_get_coor_using_mask_returns_policy_dtype(self):
        '''
        get_coor_using_mask currently returns coordinates as single precision
        '''
        #
        molecule = system.Molecule(0)
        molecule.read_pdb(PDB_DATA_PATH + '1ATM.pdb')
        mask = numpy.ones(molecule.natoms(), int)
        #
        error, coor = molecule.get_coor_using_mask(0, mask)
        #
        self.assertEqual(error, [])
        self.assertEqual(coor.dtype, config.COORD_DTYPE)

    def test_copy_molecule_using_mask_coordinates_are_policy_dtype(self):
        '''
        copy_molecule_using_mask currently stores coordinates as single precision
        '''
        #
        molecule = system.Molecule(0)
        copied = system.Molecule(1)
        molecule.read_pdb(PDB_DATA_PATH + '1ATM.pdb')
        mask = numpy.ones(molecule.natoms(), int)
        #
        error = molecule.copy_molecule_using_mask(copied, mask, 0)
        #
        self.assertEqual(error, [])
        self.assertEqual(copied.coor().dtype, config.COORD_DTYPE)

    def test_merge_two_molecules_coordinates_are_policy_dtype(self):
        '''
        merge_two_molecules stores coordinates as single precision
        '''
        #
        molecule_1 = system.Molecule_Maker(1)
        molecule_2 = system.Molecule_Maker(1)
        merged = system.Molecule(2)
        #
        error = merged.merge_two_molecules(molecule_1, molecule_2)
        #
        self.assertEqual(error, [])
        self.assertEqual(merged.coor().dtype, config.COORD_DTYPE)


if __name__ == '__main__':
    unittest.main()
