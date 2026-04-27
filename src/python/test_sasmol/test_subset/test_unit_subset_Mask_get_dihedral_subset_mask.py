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


class Test_unit_subset_Mask_get_dihedral_subset_mask(unittest.TestCase):

    def setUp(self):
        self.o = system.Molecule(0)

    def assert_list_almost_equal(self, a, b, places=5):
        if len(a) != len(b):
            raise TypeError
        for i in range(len(a)):
            if isinstance(a[i], (int, float, numpy.generic)):
                if numpy.isnan(a[i]) and numpy.isnan(b[i]):
                    continue
                self.assertAlmostEqual(a[i], b[i], places)
            else:
                self.assert_list_almost_equal(a[i], b[i], places)

    def configure_molecule(self, names, resids):
        self.o.setNatoms(len(names))
        self.o.setName(names)
        self.o.setResid(numpy.array(resids, numpy.longlong))

    def protein_expected_mask(self, flexible_residues, resids, names):
        expected = []
        for residue in flexible_residues:
            row = []
            for atom_resid, atom_name in zip(resids, names):
                include = (
                    (atom_resid == residue - 1 and atom_name == 'C') or
                    (atom_resid == residue and atom_name in ['N', 'CA', 'C']) or
                    (atom_resid == residue + 1 and atom_name == 'N')
                )
                row.append(1 if include else 0)
            expected.append(row)
        return expected

    def rna_expected_mask(self, flexible_residues, resids, names):
        expected = []
        for residue in flexible_residues:
            row = []
            for atom_resid, atom_name in zip(resids, names):
                include = (
                    (atom_resid == residue - 1 and atom_name == "O3'") or
                    (atom_resid == residue and atom_name in
                     ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]) or
                    (atom_resid == residue + 1 and atom_name == 'P') or
                    (atom_resid == residue + 1 and atom_name == "O5'")
                )
                row.append(1 if include else 0)
            expected.append(row)
        return expected

    def run_case(self, names, resids, flexible_residues, molecule_type, expected):
        self.configure_molecule(names, resids)
        result = self.o.get_dihedral_subset_mask(
            flexible_residues, molecule_type)
        self.assertTrue(isinstance(result, numpy.ndarray))
        self.assert_list_almost_equal(result.tolist(), expected)

    def test_single_residue_nomask(self):
        names = ['N', 'CA', 'C', 'O', 'CB']
        resids = [1] * len(names)
        self.run_case(names, resids, [], 0, [])

    def test_single_residue_mask(self):
        names = ['N', 'CA', 'C', 'O', 'CB']
        resids = [1] * len(names)
        expected = self.protein_expected_mask([1], resids, names)
        self.run_case(names, resids, [1], 0, expected)

    def test_three_residues_nomask(self):
        names = ['N', 'CA', 'C', 'O', 'CB', 'CG'] * 3
        resids = ([1] * 6) + ([2] * 6) + ([3] * 6)
        self.run_case(names, resids, [], 0, [])

    def test_three_residues_mask_first(self):
        names = ['N', 'CA', 'C', 'O', 'CB', 'CG'] * 3
        resids = ([1] * 6) + ([2] * 6) + ([3] * 6)
        expected = self.protein_expected_mask([1], resids, names)
        self.run_case(names, resids, [1], 0, expected)

    def test_three_residues_mask_second(self):
        names = ['N', 'CA', 'C', 'O', 'CB', 'CG'] * 3
        resids = ([1] * 6) + ([2] * 6) + ([3] * 6)
        expected = self.protein_expected_mask([2], resids, names)
        self.run_case(names, resids, [2], 0, expected)

    def test_three_residues_mask_all(self):
        names = ['N', 'CA', 'C', 'O', 'CB', 'CG'] * 3
        resids = ([1] * 6) + ([2] * 6) + ([3] * 6)
        flexible_residues = [1, 2, 3]
        expected = self.protein_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 0, expected)

    def test_500_residues_nomask(self):
        residue_names = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']
        names = residue_names * 500
        resids = []
        for residue in range(500):
            resids.extend([residue] * len(residue_names))
        self.run_case(names, resids, [], 0, [])

    def test_500_residues_mask_number_100(self):
        residue_names = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']
        names = residue_names * 500
        resids = []
        for residue in range(500):
            resids.extend([residue] * len(residue_names))
        flexible_residues = [100]
        expected = self.protein_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 0, expected)

    def test_500_residues_mask_random_10(self):
        residue_names = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']
        names = residue_names * 500
        resids = []
        for residue in range(500):
            resids.extend([residue] * len(residue_names))
        flexible_residues = [123, 12, 90, 399, 1, 89, 221, 78, 91, 129]
        expected = self.protein_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 0, expected)

    def test_500_residues_mask_all(self):
        residue_names = ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']
        names = residue_names * 500
        resids = []
        for residue in range(500):
            resids.extend([residue] * len(residue_names))
        flexible_residues = list(range(500))
        expected = self.protein_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 0, expected)

    def test_1_residue_rna_nomask(self):
        names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
        resids = [1] * len(names)
        self.run_case(names, resids, [], 1, [])

    def test_1_residue_rna_mask(self):
        names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
        resids = [1] * len(names)
        expected = self.rna_expected_mask([1], resids, names)
        self.run_case(names, resids, [1], 1, expected)

    def test_5_residues_rna_mask_random(self):
        residue_names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
        names = residue_names * 5
        resids = []
        for residue in range(1, 6):
            resids.extend([residue] * len(residue_names))
        flexible_residues = [2, 3]
        expected = self.rna_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 1, expected)

    def test_5_residues_rna_mask_all(self):
        residue_names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
        names = residue_names * 5
        resids = []
        for residue in range(1, 6):
            resids.extend([residue] * len(residue_names))
        flexible_residues = [1, 2, 3, 4, 5]
        expected = self.rna_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 1, expected)

    def test_50_residues_rna_mask_all(self):
        residue_names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]
        names = residue_names * 50
        resids = []
        for residue in range(1, 51):
            resids.extend([residue] * len(residue_names))
        flexible_residues = list(range(1, 51))
        expected = self.rna_expected_mask(flexible_residues, resids, names)
        self.run_case(names, resids, flexible_residues, 1, expected)


if __name__ == '__main__':
    main()
