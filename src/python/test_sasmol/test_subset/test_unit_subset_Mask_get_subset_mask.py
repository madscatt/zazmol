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

import numpy

import sasmol.system as system


class Test_unit_subset_Mask_get_subset_mask(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def set_parameters(self, index, name, loc, resname, chain, resid, rescode,
                       occupancy, beta, segname, element, charge, moltype,
                       residue_flag):
        natoms = len(index)
        self.o.setIndex(index)
        self.o.setName(name)
        self.o.setLoc(loc)
        self.o.setResname(resname)
        self.o.setChain(chain)
        self.o.setResid(resid)
        self.o.setRescode(rescode)
        self.o.setOccupancy(occupancy)
        self.o.setBeta(beta)
        self.o.setSegname(segname)
        self.o.setElement(element)
        self.o.setCharge(charge)
        self.o.setMoltype(moltype)
        self.o.setResidue_flag(residue_flag)
        self.o.setNatoms(natoms)

    def configure_simple_molecule(self, names, resids, elements=None,
                                  moltypes=None):
        natoms = len(names)
        if elements is None:
            elements = ['C'] * natoms
        if moltypes is None:
            moltypes = ['protein'] * natoms
        self.set_parameters(
            index=list(range(1, natoms + 1)),
            name=names,
            loc=[' '] * natoms,
            resname=['ALA'] * natoms,
            chain=['A'] * natoms,
            resid=resids,
            rescode=[' '] * natoms,
            occupancy=[1.0] * natoms,
            beta=[10.0] * natoms,
            segname=['SEG'] * natoms,
            element=elements,
            charge=[''] * natoms,
            moltype=moltypes,
            residue_flag=[0] * natoms)

    def test_null(self):
        self.configure_simple_molecule([], [])
        error, mask = self.o.get_subset_mask('')
        self.assertEqual(len(error) > 0, True)
        self.assertEqual(mask, [])

    def test_single_unmasked_due_to_no_match(self):
        self.configure_simple_molecule(['CA'], [1])
        error, mask = self.o.get_subset_mask('resid[i] == 3')
        self.assertEqual(len(error) > 0, True)
        self.assertEqual(mask, [0])

    def test_single_unmasked_due_to_syntax_error(self):
        self.configure_simple_molecule(['CA'], [1])
        error, mask = self.o.get_subset_mask(
            'name[i] == "CA" and resid[i] = 2')
        self.assertEqual(len(error) > 0, True)
        self.assertEqual(mask, [])

    def test_single_unmasked_due_to_name_error(self):
        self.configure_simple_molecule(['CA'], [1])
        error, mask = self.o.get_subset_mask('resid[js] == 3')
        self.assertEqual(len(error) > 0, True)
        self.assertEqual(mask, [])

    def test_single_masked_by_name(self):
        self.configure_simple_molecule(['CA'], [1])
        error, mask = self.o.get_subset_mask('name[i] == "CA"')
        self.assertEqual(error, [])
        self.assertTrue(isinstance(mask, numpy.ndarray))
        self.assertEqual(list(mask), [1])

    def test_single_masked_by_name_and_resid(self):
        self.configure_simple_molecule(['CA'], [1])
        error, mask = self.o.get_subset_mask(
            'name[i] == "CA" and resid[i] == 1')
        self.assertEqual(error, [])
        self.assertEqual(list(mask), [1])

    def test_three_atoms_none_masked(self):
        self.configure_simple_molecule(['CA', 'N', 'C'], [1, 1, 1],
                                       elements=['C', 'N', 'C'])
        error, mask = self.o.get_subset_mask('name[i] == "B"')
        self.assertEqual(len(error) > 0, True)
        self.assertEqual(mask, [0, 0, 0])

    def test_three_atoms_one_masked(self):
        self.configure_simple_molecule(['CA', 'N', 'C'], [1, 1, 1],
                                       elements=['C', 'N', 'C'])
        error, mask = self.o.get_subset_mask('name[i] == "CA"')
        self.assertEqual(error, [])
        self.assertEqual(list(mask), [1, 0, 0])

    def test_three_atoms_all_masked(self):
        self.configure_simple_molecule(['CA', 'N', 'C'], [1, 1, 1],
                                       elements=['C', 'N', 'C'])
        error, mask = self.o.get_subset_mask('resid[i] == 1')
        self.assertEqual(error, [])
        self.assertEqual(list(mask), [1, 1, 1])

    def test_filter_uses_element_and_moltype(self):
        self.configure_simple_molecule(
            ['CA', 'P', 'C'],
            [1, 2, 3],
            elements=['C', 'P', 'C'],
            moltypes=['protein', 'rna', 'protein'])
        error, mask = self.o.get_subset_mask(
            'element[i] == "P" and moltype[i] == "rna"')
        self.assertEqual(error, [])
        self.assertEqual(list(mask), [0, 1, 0])

    def test_named_basis_filter_all(self):
        error, basis_filter = self.o.named_basis_filter('all')
        self.assertEqual(error, [])
        self.assertEqual(basis_filter, 'not name[i] == None')

    def test_named_subset_mask_all(self):
        self.configure_simple_molecule(['H1', 'CA', 'N'], [1, 1, 1])

        error, mask = self.o.get_named_subset_mask('all')

        self.assertEqual(error, [])
        self.assertEqual(list(mask), [1, 1, 1])

    def test_named_subset_mask_heavy(self):
        self.configure_simple_molecule(['H1', 'CA', 'N', 'HA'], [1, 1, 1, 1])

        error, mask = self.o.get_named_subset_mask('heavy')

        self.assertEqual(error, [])
        self.assertEqual(list(mask), [0, 1, 1, 0])

    def test_named_subset_mask_backbone_is_not_supported(self):
        self.configure_simple_molecule(['N', 'CA', 'C', 'O'], [1, 1, 1, 1])

        error, mask = self.o.get_named_subset_mask('backbone')

        self.assertEqual(len(error) > 0, True)
        self.assertEqual(mask, [])

    def test_named_basis_filter_requires_string(self):
        error, basis_filter = self.o.named_basis_filter(None)
        self.assertEqual(len(error) > 0, True)
        self.assertEqual(basis_filter, '')


if __name__ == '__main__':
    main()
