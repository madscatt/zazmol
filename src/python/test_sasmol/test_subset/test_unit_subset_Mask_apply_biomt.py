"""Regression tests for explicit, selected BIOMT coordinate edits."""
import copy
import unittest
from unittest import mock

import numpy as np
import sasmol.config as config
import sasmol.system as system


def molecule():
    """Complete molecule with distinct frames and noncontiguous chain groups."""
    mol = system.Molecule()
    fields = dict(Atom=['ATOM'] * 4, Index=np.arange(1, 5),
                  Original_index=np.arange(1, 5), Name=['N', 'CA', 'C', 'O'],
                  Loc=[' '] * 4, Resname=['ALA'] * 4,
                  Chain=['A', 'B', 'A', 'B'], Resid=np.arange(1, 5),
                  Original_resid=np.arange(1, 5), Rescode=[' '] * 4,
                  Occupancy=['1.00'] * 4, Beta=['0.00'] * 4,
                  Segname=['S1', 'S2', 'S1', 'S2'], Element=['N', 'C', 'C', 'O'],
                  Charge=[' '] * 4, Moltype=['protein'] * 4,
                  Residue_flag=[False] * 4, Charmm_type=[],
                  Atom_charge=np.zeros(4), Atom_vdw=np.zeros(4), Conect={})
    for field, value in fields.items():
        getattr(mol, 'set' + field)(value)
    mol.setNatoms(4)
    mol.setNumber_of_frames(2)
    first = np.array([[1, 2, 3], [-2, 4, 5], [3, -1, 2], [4, 3, -2]],
                     dtype=config.COORD_DTYPE)
    mol.setCoor(np.stack([first, first + 10]))
    return mol


class TestApplyBiomt(unittest.TestCase):
    def test_float64_destination_preserves_small_translation(self):
        # Extracted coordinates use COORD_DTYPE, but transform results must
        # retain the precision of an existing CALC_DTYPE destination.
        translation = np.array([1e-8, 2e-8, 3e-8], dtype=config.CALC_DTYPE)
        for frame in (0, 1):
            for indices in ([0, 1, 2, 3], [1, 3]):
                with self.subTest(frame=frame, indices=indices):
                    mol = molecule()
                    mol.setCoor(mol.coor().astype(config.CALC_DTYPE))
                    coordinates = mol.coor()
                    expected = coordinates.copy()
                    expected[frame, indices] += translation
                    mol.apply_biomt(frame, 'i in ' + repr(indices),
                                    np.eye(3, dtype=config.CALC_DTYPE), translation)
                    self.assertIs(mol.coor(), coordinates)
                    self.assertEqual(mol.coor().dtype, np.dtype(config.CALC_DTYPE))
                    np.testing.assert_array_equal(mol.coor(), expected)

    def test_selected_indices_and_frames(self):
        for frame in (0, 1, -1, -2, np.int64(1)):
            for indices in ([1], [1, 3], [0, 1, 2, 3]):
                with self.subTest(frame=frame, indices=indices):
                    mol = molecule()
                    original = mol.coor()
                    expected = original.copy()
                    expected[frame, indices] += [10, 20, 30]
                    descriptors = copy.deepcopy(mol.__dict__)
                    result = mol.apply_biomt(frame, 'i in ' + repr(indices),
                                             np.eye(3), [10, 20, 30], future=True)
                    self.assertIsNone(result)
                    self.assertIs(mol.coor(), original)
                    np.testing.assert_array_equal(mol.coor(), expected)
                    self.assertEqual(mol.coor().dtype, expected.dtype)
                    for key, value in descriptors.items():
                        if key != '_coor':
                            np.testing.assert_equal(mol.__dict__[key], value)

    def test_rotation_translation_hand_calculated(self):
        mol = molecule()
        before = mol.coor().copy()
        rotation = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        mol.apply_biomt(0, 'i < 2', rotation, [10, 20, 30])
        before[0, :2] = [[8, 21, 33], [6, 18, 35]]
        np.testing.assert_array_equal(mol.coor(), before)

    def test_identity_and_nonorthogonal_coefficients(self):
        mol = molecule()
        original = mol.coor().copy()
        mol.apply_biomt(0, 'True', np.eye(3), np.zeros(3))
        np.testing.assert_array_equal(mol.coor(), original)
        mol.apply_biomt(0, 'i == 0', np.diag([2., 3., 4.]), [1, 2, 3])
        original[0, 0] = [3, 8, 15]
        np.testing.assert_array_equal(mol.coor(), original)

    def test_translation_representations(self):
        for count in (1, 3, 4):
            for shape in ((3,), (1, 3), (3, 1)):
                with self.subTest(count=count, shape=shape):
                    mol = molecule()
                    expected = mol.coor().copy()
                    expected[0, :count] += [1, 2, 3]
                    mol.apply_biomt(0, 'i < ' + str(count), np.eye(3),
                                    np.array([1., 2., 3.]).reshape(shape))
                    np.testing.assert_array_equal(mol.coor(), expected)

    def test_selection_errors_are_reported_before_extraction(self):
        for selection in ('False', 'bad[', 'missing[i]', 'i == 0 or missing[i]'):
            with self.subTest(selection=selection):
                mol = molecule()
                before = mol.coor().copy()
                with mock.patch.object(mol, 'get_coor_using_mask') as extract:
                    with self.assertRaisesRegex(ValueError, 'selection|filter'):
                        mol.apply_biomt(0, selection, np.eye(3), np.zeros(3))
                    extract.assert_not_called()
                np.testing.assert_array_equal(mol.coor(), before)

    def test_invalid_frames_preserve_coordinates(self):
        for frame in (2, -3):
            mol = molecule()
            before = mol.coor().copy()
            with self.assertRaises(IndexError):
                mol.apply_biomt(frame, 'True', np.eye(3), np.zeros(3))
            np.testing.assert_array_equal(mol.coor(), before)

    def test_calculation_failures_preserve_coordinates(self):
        for rotation, translation in ((np.zeros((2, 2)), np.zeros(3)),
                                      (np.eye(3), np.zeros(2))):
            mol = molecule()
            before = mol.coor().copy()
            with self.assertRaises(ValueError):
                mol.apply_biomt(0, 'True', rotation, translation)
            np.testing.assert_array_equal(mol.coor(), before)

    def test_extraction_errors_and_exceptions(self):
        mol = molecule()
        before = mol.coor().copy()
        with mock.patch.object(mol, 'get_coor_using_mask',
                               return_value=(['extraction diagnostic'], None)):
            with self.assertRaisesRegex(ValueError, 'extraction diagnostic'):
                mol.apply_biomt(0, 'True', np.eye(3), np.zeros(3))
        with mock.patch.object(mol, 'get_coor_using_mask',
                               side_effect=RuntimeError('extraction exception')):
            with self.assertRaisesRegex(RuntimeError, 'extraction exception'):
                mol.apply_biomt(0, 'True', np.eye(3), np.zeros(3))
        np.testing.assert_array_equal(mol.coor(), before)


if __name__ == '__main__':
    unittest.main()
