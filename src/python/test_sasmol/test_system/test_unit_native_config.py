import os
import shutil
import tempfile
import unittest

import sasmol.native_config as native_config


class Test_unit_native_config(unittest.TestCase):

    def setUp(self):
        self.original_file = native_config.__file__
        self.temp_dir = tempfile.mkdtemp()
        package_dir = os.path.join(self.temp_dir, "sasmol")
        native_dir = os.path.join(package_dir, "native")
        os.makedirs(os.path.join(native_dir, "include", "sasmol"))
        os.makedirs(os.path.join(native_dir, "lib"))
        open(os.path.join(package_dir, "native_config.py"), "w").close()
        open(os.path.join(native_dir, "include", "sasmol", "file_io.hpp"), "w").close()
        open(os.path.join(native_dir, "lib", "libsasmol.a"), "w").close()
        native_config.__file__ = os.path.join(package_dir, "native_config.py")

    def tearDown(self):
        native_config.__file__ = self.original_file
        shutil.rmtree(self.temp_dir)

    def test_paths_point_to_packaged_native_payload(self):
        include_dirs = native_config.include_dirs()
        library_dirs = native_config.library_dirs()
        extra_objects = native_config.extra_objects()

        self.assertEqual(1, len(include_dirs))
        self.assertEqual(1, len(library_dirs))
        self.assertEqual(1, len(extra_objects))
        self.assertTrue(include_dirs[0].endswith(os.path.join("native", "include")))
        self.assertTrue(library_dirs[0].endswith(os.path.join("native", "lib")))
        self.assertTrue(extra_objects[0].endswith(
            os.path.join("native", "lib", "libsasmol.a")))

    def test_modern_dcd_capability_is_advertised_without_legacy_handle_abi(self):
        self.assertGreaterEqual(native_config.native_api_version(), 1)
        self.assertTrue(native_config.has_modern_dcd_reader_api())
        self.assertTrue(native_config.has_capability("dcd_read_next_frame_coordinates"))
        self.assertFalse(native_config.has_dcd_read_handle_api())


if __name__ == "__main__":
    unittest.main()
