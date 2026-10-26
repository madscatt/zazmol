"""Real copy/transform compatibility and error propagation regressions."""
import copy
import unittest
from unittest import mock

import numpy as np
import sasmol.config as config
import sasmol.system as system
from sasmol.test_sasmol.test_subset.test_unit_subset_Mask_apply_biomt import molecule


class TestCopyApplyBiomt(unittest.TestCase):
    def test_copy_retains_coordinate_storage_dtype_from_float64_source(self):
        source, destination = molecule(), system.Molecule()
        source.setCoor(source.coor().astype(config.CALC_DTYPE))
        before = source.coor().copy()
        source.copy_apply_biomt(destination, 1, 'i == 1 or i == 3',
                                np.eye(3), [1e-8, 2e-8, 3e-8])
        self.assertEqual(destination.coor().dtype, np.dtype(config.COORD_DTYPE))
        expected = (before[1, [1, 3]] + [1e-8, 2e-8, 3e-8]).astype(config.COORD_DTYPE)
        np.testing.assert_array_equal(destination.coor(), expected[None, :, :])
        np.testing.assert_array_equal(source.coor(), before)

    def test_source_selection_and_source_frame(self):
        for frame in (0, 1, -1):
            for selection, indices in (('i == 1 or i == 3', [1, 3]),
                                       ('i == 3', [3]),
                                       ("chain[i] == 'B'", [1, 3]),
                                       ('True', [0, 1, 2, 3])):
                with self.subTest(frame=frame, selection=selection):
                    source, destination = molecule(), system.Molecule()
                    before = copy.deepcopy(source.__dict__)
                    with mock.patch.object(source, 'get_subset_mask',
                                           wraps=source.get_subset_mask) as select:
                        result = source.copy_apply_biomt(
                            destination, frame, selection, np.eye(3), [1, 2, 3],
                            future=True)
                        self.assertEqual(select.call_count, 1)
                    self.assertIsNone(result)
                    np.testing.assert_array_equal(
                        destination.coor(),
                        before['_coor'][frame, indices][None, :, :] + [1, 2, 3])
                    self.assertEqual(destination.natoms(), len(indices))
                    self.assertEqual(destination.coor().dtype, source.coor().dtype)
                    self.assertEqual(destination.name(), [source.name()[i] for i in indices])
                    self.assertEqual(destination.chain(), [source.chain()[i] for i in indices])
                    for key, value in before.items():
                        np.testing.assert_equal(source.__dict__[key], value)

    def test_copy_rotation_and_translation_representations(self):
        for shape in ((3,), (1, 3), (3, 1)):
            source, destination = molecule(), system.Molecule()
            source.copy_apply_biomt(destination, 0, 'i < 2',
                                    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                                    np.array([10., 20., 30.]).reshape(shape))
            np.testing.assert_array_equal(destination.coor(),
                                          [[[8, 21, 33], [6, 18, 35]]])

    def test_existing_caller_pattern(self):
        for selection in ("chain[i] == 'B'", "segname[i] == 'S2'"):
            source, destination = molecule(), system.Molecule()
            before = source.coor().copy()
            errors, mask = source.get_subset_mask(selection)
            self.assertEqual(errors, [])
            self.assertEqual(source.copy_molecule_using_mask(destination, mask, 1), [])
            destination.apply_biomt(0, selection,
                                    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
                                    [10, 20, 30])
            np.testing.assert_array_equal(destination.coor(),
                                          [[[-4, 28, 45], [-3, 34, 38]]])
            np.testing.assert_array_equal(source.coor(), before)
            self.assertEqual(destination.name(), ['CA', 'O'])

    def test_selection_errors_do_not_copy(self):
        for selection in ('False', 'bad[', 'i == 0 or missing[i]'):
            source, destination = molecule(), molecule()
            before = destination.coor().copy()
            with mock.patch.object(source, 'copy_molecule_using_mask') as copier:
                with self.assertRaisesRegex(ValueError, 'selection|filter'):
                    source.copy_apply_biomt(destination, 0, selection, np.eye(3), [0, 0, 0])
                copier.assert_not_called()
            np.testing.assert_array_equal(destination.coor(), before)

    def test_bad_operation_and_frame_do_not_copy(self):
        for frame, rotation, translation, exception in (
                (2, np.eye(3), [0, 0, 0], IndexError),
                (0, np.zeros((2, 2)), [0, 0, 0], ValueError)):
            source, destination = molecule(), molecule()
            before = destination.coor().copy()
            with mock.patch.object(source, 'copy_molecule_using_mask') as copier:
                with self.assertRaises(exception):
                    source.copy_apply_biomt(destination, frame, 'True', rotation, translation)
                copier.assert_not_called()
            np.testing.assert_array_equal(destination.coor(), before)

    def test_copy_failure_does_not_transform_partial_destination(self):
        for raises in (False, True):
            source, destination = molecule(), molecule()
            before = copy.deepcopy(source.__dict__)
            def fail_copy(other, mask, frame):
                other.setCoor(np.array([[[9., 8., 7.]]], dtype=source.coor().dtype))
                other.setName(['PARTIAL'])
                if raises:
                    raise RuntimeError('copy exception')
                return ['copy diagnostic']
            with mock.patch.object(source, 'copy_molecule_using_mask', side_effect=fail_copy):
                with mock.patch.object(destination, 'apply_biomt') as transform:
                    with self.assertRaisesRegex(RuntimeError if raises else ValueError,
                                                'copy exception' if raises else 'copy diagnostic'):
                        source.copy_apply_biomt(destination, 0, 'True', np.eye(3), [1, 2, 3])
                    transform.assert_not_called()
            # The copied failure state is deliberately not restored; no transform follows it.
            np.testing.assert_array_equal(destination.coor(), [[[9., 8., 7.]]])
            for key, value in before.items():
                np.testing.assert_equal(source.__dict__[key], value)


if __name__ == '__main__':
    unittest.main()
