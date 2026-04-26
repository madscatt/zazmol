'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
'''

from unittest import main
import unittest

import numpy

import sasmol.mask as mask


class Test_unit_subset_Mask_get_mask_array_extension(unittest.TestCase):

    def protein_case(self):
        names = ['N', 'CA', 'C', 'O'] * 3
        resids = [1] * 4 + [2] * 4 + [3] * 4
        flexible_residues = [2]
        expected = [[
            0, 0, 1, 0,
            1, 1, 1, 0,
            1, 0, 0, 0,
        ]]
        return names, resids, flexible_residues, expected

    def rna_case(self):
        residue_names = ["P", "O5'", "C5'", "C4'", "C3'", "O3'", "C1'"]
        names = residue_names * 3
        resids = [1] * 7 + [2] * 7 + [3] * 7
        flexible_residues = [2]
        expected = [[
            0, 0, 0, 0, 0, 1, 0,
            1, 1, 1, 1, 1, 1, 0,
            1, 1, 0, 0, 0, 0, 0,
        ]]
        return names, resids, flexible_residues, expected

    def call_mask(self, names, resids, flexible_residues, molecule_type,
                  farray_dtype=numpy.longlong, resid_dtype=numpy.longlong,
                  flexible_dtype=numpy.longlong):
        farray = numpy.zeros(
            (len(flexible_residues), len(names)), dtype=farray_dtype)
        resid = numpy.array(resids, dtype=resid_dtype)
        flexible = numpy.array(flexible_residues, dtype=flexible_dtype)

        result = mask.get_mask_array(
            farray, names, resid, flexible, len(set(resids)), molecule_type)

        return result, farray

    def test_protein_mask_mutates_output_array_in_place(self):
        names, resids, flexible_residues, expected = self.protein_case()

        result, farray = self.call_mask(names, resids, flexible_residues, 0)

        self.assertIsNone(result)
        self.assertEqual(farray.dtype, numpy.dtype(numpy.longlong))
        self.assertEqual(farray.tolist(), expected)

    def test_rna_mask_mutates_output_array_in_place(self):
        names, resids, flexible_residues, expected = self.rna_case()

        result, farray = self.call_mask(names, resids, flexible_residues, 1)

        self.assertIsNone(result)
        self.assertEqual(farray.dtype, numpy.dtype(numpy.longlong))
        self.assertEqual(farray.tolist(), expected)

    def test_multiple_flexible_residues_create_multiple_rows(self):
        names = ['N', 'CA', 'C', 'O'] * 4
        resids = [1] * 4 + [2] * 4 + [3] * 4 + [4] * 4
        flexible_residues = [2, 3]
        expected = [
            [0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
        ]

        result, farray = self.call_mask(names, resids, flexible_residues, 0)

        self.assertIsNone(result)
        self.assertEqual(farray.tolist(), expected)

    def test_empty_flexible_residues_leaves_empty_output(self):
        names, resids, _, _ = self.protein_case()

        result, farray = self.call_mask(names, resids, [], 0)

        self.assertIsNone(result)
        self.assertEqual(farray.shape, (0, len(names)))

    def test_invalid_molecule_type_leaves_output_zeroed(self):
        names, resids, flexible_residues, _ = self.protein_case()

        result, farray = self.call_mask(names, resids, flexible_residues, 99)

        self.assertIsNone(result)
        self.assertEqual(farray.tolist(), [[0] * len(names)])

    def test_name_argument_must_be_a_list(self):
        names, resids, flexible_residues, _ = self.protein_case()
        farray = numpy.zeros(
            (len(flexible_residues), len(names)), dtype=numpy.longlong)

        with self.assertRaisesRegex(TypeError, 'not a list'):
            mask.get_mask_array(
                farray, tuple(names), numpy.array(resids, numpy.longlong),
                numpy.array(flexible_residues, numpy.longlong),
                len(set(resids)), 0)

    def test_name_list_members_must_be_strings(self):
        names, resids, flexible_residues, _ = self.protein_case()
        names[1] = 42
        farray = numpy.zeros(
            (len(flexible_residues), len(names)), dtype=numpy.longlong)

        with self.assertRaisesRegex(TypeError, 'list must contain strings'):
            mask.get_mask_array(
                farray, names, numpy.array(resids, numpy.longlong),
                numpy.array(flexible_residues, numpy.longlong),
                len(set(resids)), 0)

    def test_output_array_must_be_two_dimensional(self):
        names, resids, flexible_residues, _ = self.protein_case()
        farray = numpy.zeros(len(names), dtype=numpy.longlong)

        with self.assertRaisesRegex(TypeError, 'Array must have 2 dimensions'):
            mask.get_mask_array(
                farray, names, numpy.array(resids, numpy.longlong),
                numpy.array(flexible_residues, numpy.longlong),
                len(set(resids)), 0)

    def test_output_array_must_be_longlong(self):
        names, resids, flexible_residues, _ = self.protein_case()
        farray = numpy.zeros(
            (len(flexible_residues), len(names)), dtype=numpy.int32)

        with self.assertRaisesRegex(TypeError, "Array of type 'long long'"):
            mask.get_mask_array(
                farray, names, numpy.array(resids, numpy.longlong),
                numpy.array(flexible_residues, numpy.longlong),
                len(set(resids)), 0)

    def test_input_residue_arrays_accept_int32(self):
        names, resids, flexible_residues, expected = self.protein_case()

        result, farray = self.call_mask(
            names, resids, flexible_residues, 0,
            resid_dtype=numpy.int32, flexible_dtype=numpy.int32)

        self.assertIsNone(result)
        self.assertEqual(farray.tolist(), expected)


if __name__ == '__main__':
    main()
