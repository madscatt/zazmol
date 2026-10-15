import pathlib
import subprocess
import sys

REPOSITORY = pathlib.Path(__file__).resolve().parents[4]
COMPONENT = "zazmol"
OWNER_VERSION = "2.0"


def load_build_support():
    import importlib.util
    spec = importlib.util.spec_from_file_location("owner_build_support", REPOSITORY / "build_support.py")
    support = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(support)
    return support


def test_build_stamp_tracks_built_source_and_marks_development_changes(tmp_path):
    import json
    support = load_build_support()
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    source = tmp_path / "source.txt"
    source.write_text("one")
    subprocess.run(["git", "-C", str(tmp_path), "add", "source.txt"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=Fixture", "-c", "user.email=fixture@example.test", "commit", "-qm", "first"], check=True)
    first = support.build_info(tmp_path, COMPONENT, OWNER_VERSION)
    assert not first["modified"]
    assert len(first["revision"]) == 40
    source.write_text("two")
    assert support.build_info(tmp_path, COMPONENT, OWNER_VERSION)["modified"]
    receipt = tmp_path / "installed/_build_info.json"
    support.write_info(receipt, first)
    assert json.loads(receipt.read_text()) == first  # Later source edits cannot change an installed record.


def test_source_archive_receipt_survives_no_git_and_unknown_stays_unknown(tmp_path):
    import json
    support = load_build_support()
    assert support.build_info(tmp_path, COMPONENT, OWNER_VERSION)["revision"] is None
    receipt = dict(schema_version=1, component_id=COMPONENT, version=OWNER_VERSION,
                   revision="b" * 40, modified=False, built_at="2026-10-15T00:00:00Z")
    (tmp_path / ".build-info.json").write_text(json.dumps(receipt))
    assert support.build_info(tmp_path, COMPONENT, OWNER_VERSION)["revision"] == "b" * 40
    receipt["component_id"] = "wrong_owner"
    (tmp_path / ".build-info.json").write_text(json.dumps(receipt))
    import pytest
    with pytest.raises(ValueError): support.build_info(tmp_path, COMPONENT, OWNER_VERSION)


def test_real_egg_and_source_archive_include_stamp_and_failed_build_has_no_stamp(tmp_path):
    import json
    import os
    import shutil
    import tarfile
    import zipfile
    shutil.copy2(REPOSITORY / "build_support.py", tmp_path / "build_support.py")
    (tmp_path / "fixturepkg").mkdir()
    (tmp_path / "fixturepkg/__init__.py").write_text("")
    (tmp_path / "README").write_text("A packaging fixture, not a scientific package.\n")
    (tmp_path / "setup.py").write_text(f'''
from setuptools import setup
from distutils.command.build import build as Build
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist
from build_support import build_commands
import os
class CheckedBuild(Build):
    def run(self):
        super().run()
        if os.environ.get("FAIL_BUILD"):
            raise RuntimeError("intentional fixture build failure")
setup(name="stamp-fixture", version="{OWNER_VERSION}", packages=["fixturepkg"],
      cmdclass=build_commands(os.getcwd(), "{COMPONENT}", "{OWNER_VERSION}",
                              "fixturepkg", CheckedBuild, build_py, sdist))
''')
    command = [sys.executable, "setup.py"]
    subprocess.run(command + ["bdist_egg", "sdist"], cwd=tmp_path, check=True, capture_output=True)
    egg = next((tmp_path / "dist").glob("*.egg"))
    with zipfile.ZipFile(egg) as archive:
        stamp = json.loads(archive.read("fixturepkg/_build_info.json"))
    assert stamp["component_id"] == COMPONENT
    assert stamp["version"] == OWNER_VERSION
    with tarfile.open(next((tmp_path / "dist").glob("*.tar.gz"))) as archive:
        assert any(n.endswith("/build_support.py") for n in archive.getnames())
        assert any(n.endswith("/.build-info.json") for n in archive.getnames())
    env = dict(os.environ, FAIL_BUILD="1")
    failed = subprocess.run(command + ["build"], cwd=tmp_path, env=env, capture_output=True)
    assert failed.returncode != 0
    assert not list((tmp_path / "build").glob("lib*/fixturepkg/_build_info.json"))
