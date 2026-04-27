#!/usr/bin/env python
"""
Build the project and run the Python test suite from the build tree.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import unittest


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
TEST_ROOT = os.path.join(REPO_ROOT, "src", "python", "test_sasmol")

SUITES = {
    "all": TEST_ROOT,
    "system": os.path.join(TEST_ROOT, "test_system"),
    "calculate": os.path.join(TEST_ROOT, "test_calculate"),
    "file_io": os.path.join(TEST_ROOT, "test_file_io"),
    "linear_algebra": os.path.join(TEST_ROOT, "test_linear_algebra"),
    "operate": os.path.join(TEST_ROOT, "test_operate"),
    "properties": os.path.join(TEST_ROOT, "test_properties"),
    "subset": os.path.join(TEST_ROOT, "test_subset"),
    "topology": os.path.join(TEST_ROOT, "test_topology"),
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build zazmol and run Python tests from the build tree."
    )
    parser.add_argument(
        "--suite",
        choices=sorted(SUITES.keys()),
        default="all",
        help="test slice to run (default: all)",
    )
    parser.add_argument(
        "--no-build",
        action="store_true",
        help="reuse an existing build directory instead of rebuilding first",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list the supported suites and exit",
    )
    return parser.parse_args()


def build_package():
    build_dir = os.path.join(REPO_ROOT, "build")
    if os.path.isdir(build_dir):
        shutil.rmtree(build_dir)

    subprocess.check_call([sys.executable, "setup.py", "build"], cwd=REPO_ROOT)


def find_build_lib():
    matches = glob.glob(os.path.join(REPO_ROOT, "build", "lib.*"))
    if not matches:
        raise RuntimeError(
            "build/lib.* was not created; run the build step first")
    return matches[0]


def run_suite(suite_name):
    suite_path = SUITES[suite_name]
    build_lib = find_build_lib()
    sys.path.insert(0, build_lib)

    loader = unittest.defaultTestLoader
    suite = loader.discover(suite_path, pattern="test_*.py")
    result = unittest.TextTestRunner().run(suite)
    return 0 if result.wasSuccessful() else 1


def main():
    args = parse_args()

    if args.list:
        for suite_name in sorted(SUITES.keys()):
            print(suite_name)
        return 0

    if not args.no_build:
        build_package()

    return run_suite(args.suite)


if __name__ == "__main__":
    raise SystemExit(main())
