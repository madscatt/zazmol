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

import sasmol.system as system


class Test_intg_system_Molecule_Maker(unittest.TestCase):

    def test_default_constructor_values(self):
        '''
        test default constructor values
        '''

        natoms = 3
        molecule = system.Molecule_Maker(natoms)

        self.assertEqual(molecule.natoms(), natoms)
        self.assertEqual(molecule.number_of_atoms(), natoms)
        self.assertEqual(molecule.atom(), ['ATOM', 'ATOM', 'ATOM'])
        self.assertEqual(list(molecule.index()), [1, 2, 3])
        self.assertEqual(molecule.name(), ['C', 'C', 'C'])
        self.assertEqual(molecule.loc(), [' ', ' ', ' '])
        self.assertEqual(molecule.resname(), ['DUM', 'DUM', 'DUM'])
        self.assertEqual(molecule.chain(), ['A', 'A', 'A'])
        self.assertEqual(list(molecule.resid()), [1, 2, 3])
        self.assertEqual(molecule.rescode(), [' ', ' ', ' '])
        self.assertEqual(molecule.coor().shape, (1, natoms, 3))
        self.assertEqual(molecule.coor().dtype, numpy.float32)
        self.assertEqual(molecule.occupancy(), ['0.00', '0.00', '0.00'])
        self.assertEqual(molecule.beta(), ['0.00', '0.00', '0.00'])
        self.assertEqual(molecule.charge(), [' ', ' ', ' '])
        self.assertEqual(molecule.segname(), ['DUM', 'DUM', 'DUM'])
        self.assertEqual(molecule.element(), ['C', 'C', 'C'])
        self.assertEqual(molecule.moltype(), ['protein', 'protein', 'protein'])
        self.assertEqual(molecule.residue_flag(), [False, False, False])
        self.assertEqual(list(molecule.original_index()), [1, 2, 3])
        self.assertEqual(list(molecule.original_resid()), [1, 2, 3])

    def test_validate_integrity_passes_for_constructor_output(self):
        '''
        test constructor output passes strict integrity validation
        '''

        molecule = system.Molecule_Maker(2, name=['N', 'CA'])

        result = molecule.validate_integrity()

        self.assertEqual(result['_name'], 2)
        self.assertEqual(result['_index'], 2)
        self.assertEqual(result['_resid'], 2)

    def test_validate_integrity_raises_for_setter_length_mismatch(self):
        '''
        test strict integrity validation catches later setter mismatches
        '''

        molecule = system.Molecule_Maker(2)
        molecule.setName(['N'])

        with self.assertRaisesRegex(ValueError, '_name has length 1'):
            molecule.validate_integrity()

    def test_scalar_constructor_values_are_broadcast(self):
        '''
        test scalar constructor values are assigned to every atom
        '''

        natoms = 2
        molecule = system.Molecule_Maker(
            natoms,
            atom='HETATM',
            name='Ar',
            loc='A',
            resname='ARG',
            chain='B',
            rescode='C',
            occupancy='1.00',
            beta='2.00',
            charge='+1',
            segname='ARG0',
            element='Ar',
            moltype='other',
            residue_flag=True)

        self.assertEqual(molecule.atom(), ['HETATM', 'HETATM'])
        self.assertEqual(molecule.name(), ['Ar', 'Ar'])
        self.assertEqual(molecule.loc(), ['A', 'A'])
        self.assertEqual(molecule.resname(), ['ARG', 'ARG'])
        self.assertEqual(molecule.chain(), ['B', 'B'])
        self.assertEqual(molecule.rescode(), ['C', 'C'])
        self.assertEqual(molecule.occupancy(), ['1.00', '1.00'])
        self.assertEqual(molecule.beta(), ['2.00', '2.00'])
        self.assertEqual(molecule.charge(), ['+1', '+1'])
        self.assertEqual(molecule.segname(), ['ARG0', 'ARG0'])
        self.assertEqual(molecule.element(), ['Ar', 'Ar'])
        self.assertEqual(molecule.moltype(), ['other', 'other'])
        self.assertEqual(molecule.residue_flag(), [True, True])

    def test_index_resid_and_coordinates_can_be_passed_directly(self):
        '''
        test constructor values that already support per-atom arrays
        '''

        index = numpy.array([10, 11], numpy.int32)
        resid = [20, 21]
        coor = numpy.array([[[1.0, 2.0, 3.0],
                             [4.0, 5.0, 6.0]]], numpy.float32)

        molecule = system.Molecule_Maker(
            2, index=index, resid=resid, coor=coor)

        self.assertIs(molecule.index(), index)
        self.assertIs(molecule.resid(), resid)
        self.assertIs(molecule.coor(), coor)
        self.assertIs(molecule.original_index(), index)
        self.assertIs(molecule.original_resid(), resid)

    def test_per_atom_constructor_values_are_accepted(self):
        '''
        test list and array values can be supplied per atom
        '''

        molecule = system.Molecule_Maker(
            3,
            atom=['ATOM', 'ATOM', 'HETATM'],
            name=['N', 'CA', 'C'],
            loc=[' ', 'A', ' '],
            resname=['GLY', 'ALA', 'HOH'],
            chain=['A', 'A', 'B'],
            rescode=[' ', ' ', 'X'],
            occupancy=['1.00', '0.50', '0.25'],
            beta=('0.00', '1.00', '2.00'),
            charge=[' ', '+1', '-1'],
            segname=['SEG1', 'SEG1', 'WAT'],
            element=numpy.array(['N', 'C', 'O']),
            moltype=['protein', 'protein', 'water'],
            residue_flag=[False, False, True])

        self.assertEqual(molecule.atom(), ['ATOM', 'ATOM', 'HETATM'])
        self.assertEqual(molecule.name(), ['N', 'CA', 'C'])
        self.assertEqual(molecule.loc(), [' ', 'A', ' '])
        self.assertEqual(molecule.resname(), ['GLY', 'ALA', 'HOH'])
        self.assertEqual(molecule.chain(), ['A', 'A', 'B'])
        self.assertEqual(molecule.rescode(), [' ', ' ', 'X'])
        self.assertEqual(molecule.occupancy(), ['1.00', '0.50', '0.25'])
        self.assertEqual(molecule.beta(), ['0.00', '1.00', '2.00'])
        self.assertEqual(molecule.charge(), [' ', '+1', '-1'])
        self.assertEqual(molecule.segname(), ['SEG1', 'SEG1', 'WAT'])
        self.assertEqual(molecule.element(), ['N', 'C', 'O'])
        self.assertEqual(
            molecule.moltype(), ['protein', 'protein', 'water'])
        self.assertEqual(molecule.residue_flag(), [False, False, True])

    def test_strings_are_broadcast_not_split(self):
        '''
        test strings remain scalar constructor values
        '''

        molecule = system.Molecule_Maker(3, name='CA', chain='A')

        self.assertEqual(molecule.name(), ['CA', 'CA', 'CA'])
        self.assertEqual(molecule.chain(), ['A', 'A', 'A'])

    def test_wrong_length_per_atom_values_raise_value_error(self):
        '''
        test per-atom values must match the number of atoms
        '''

        values = {
            'atom': ['ATOM'],
            'index': [1],
            'name': ['N'],
            'loc': [' '],
            'resname': ['GLY'],
            'chain': ['A'],
            'resid': [1],
            'rescode': [' '],
            'occupancy': ['1.00'],
            'beta': ['0.00'],
            'charge': [' '],
            'segname': ['SEG1'],
            'element': ['N'],
            'moltype': ['protein'],
            'residue_flag': [False],
        }

        for field_name, value in values.items():
            with self.assertRaisesRegex(ValueError, field_name):
                system.Molecule_Maker(2, **{field_name: value})


if __name__ == '__main__':
    unittest.main()
