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


def test_install_completion_is_optional_and_bound_to_its_environment(tmp_path, monkeypatch):
    import json
    import pytest
    support = load_build_support()
    monkeypatch.setattr(support.sys, "prefix", str(tmp_path))
    expected = tmp_path / "verified build.json"
    expected.write_text("{}")
    # No configuration means a standalone install, with no application needed.
    support.complete_install(COMPONENT, expected)
    publisher = tmp_path / "publisher with spaces.py"
    publisher.write_text("import json, pathlib, sys\np=pathlib.Path(__file__).parent\n"
                         "(p/'arguments.json').write_text(json.dumps(sys.argv[1:]))\n"
                         "raise SystemExit(8 if (p/'fail').exists() else 0)\n")
    config_path = tmp_path / "etc/sassie/version-stack.json"
    config_path.parent.mkdir(parents=True)
    config = dict(schema_version=1, python_prefix=str(tmp_path), publisher=str(publisher))
    config_path.write_text(json.dumps(config))
    support.complete_install(COMPONENT, expected)
    assert json.loads((tmp_path / "arguments.json").read_text()) == [
        "complete-install", "--component", COMPONENT, "--expected-stamp", str(expected), "--verified"]
    (tmp_path / "fail").touch()
    with pytest.raises(RuntimeError, match="Installed successfully; version-stack update failed"):
        support.complete_install(COMPONENT, expected)
    (tmp_path / "arguments.json").unlink()
    config["python_prefix"] = str(tmp_path / "another environment")
    config_path.write_text(json.dumps(config))
    with pytest.raises(RuntimeError, match="invalid version-stack binding"):
        support.complete_install(COMPONENT, expected)
    assert not (tmp_path / "arguments.json").exists()
    config_path.unlink()
    config_path.symlink_to(tmp_path / "missing binding")
    with pytest.raises(RuntimeError, match="Installed successfully; version-stack update failed"):
        support.complete_install(COMPONENT, expected)


def test_install_verification_and_completion_order_and_staging(tmp_path, monkeypatch):
    from types import SimpleNamespace
    import pytest
    import sysconfig
    support = load_build_support()
    calls = []
    class Installer:
        root = home = None
        user = dry_run = False
        prefix = sys.prefix
        install_lib = sysconfig.get_path("purelib")
        def run(self):
            calls.append("install")
            if failure == "install":
                raise RuntimeError("install failed")
        def get_finalized_command(self, name):
            return SimpleNamespace(build_lib=tmp_path)
    def verify(expected):
        calls.append("verify")
        assert expected == tmp_path / "sasmol/_build_info.json"
        if failure == "verify":
            raise RuntimeError("verify failed")
    def complete(component, expected):
        calls.append("complete")
        assert component == "zazmol"
        if failure == "complete":
            raise RuntimeError("completion failed")
    monkeypatch.setattr(support, "verify_installed_zazmol", verify)
    monkeypatch.setattr(support, "complete_install", complete)
    for failure in ("install", "verify", "complete", "none"):
        calls.clear()
        command = support.verified_install_command(Installer)()
        if failure == "none":
            command.run()
        else:
            with pytest.raises(RuntimeError): command.run()
        expected = ["install", "verify", "complete"]
        assert calls == (expected if failure == "none" else expected[:expected.index(failure) + 1])
    for change in ({"root": str(tmp_path)}, {"prefix": str(tmp_path)}, {"home": str(tmp_path)},
                   {"user": True}, {"dry_run": True}, {"install_lib": str(tmp_path)}):
        calls.clear()
        command = support.verified_install_command(Installer)()
        for key, value in change.items(): setattr(command, key, value)
        command.run()
        assert calls == ["install"]


def test_normal_install_updates_bound_stack_and_staging_does_not(tmp_path):
    """Exercise real setuptools install with a tiny metadata fixture, not science."""
    import json
    import shutil
    import venv
    import pytest
    app = REPOSITORY.parent / "genapp_zazzie"
    publisher = app / "tools/source_revision_metadata.py"
    if not publisher.is_file():
        pytest.skip("optional paired application checkout is unavailable")
    prefix = tmp_path / "fixture environment"
    venv.EnvBuilder(with_pip=False, system_site_packages=True).create(prefix)
    python = prefix / "bin/python"
    fixture = tmp_path / "fixture source"
    fixture.mkdir()
    shutil.copy2(REPOSITORY / "build_support.py", fixture / "build_support.py")
    (fixture / "sasmol").mkdir()
    (fixture / "sasmol/__init__.py").write_text("")
    # Native behavior is checked independently. This fixture verifies the real
    # install-command dispatch and the actual application completion writer.
    (fixture / "setup.py").write_text('''
import os, subprocess, sys
from pathlib import Path
from setuptools import setup
from distutils.command.build import build
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist
from setuptools.command.install import install
import build_support
def verify_fixture(expected):
    script = "import json,sys; from importlib import resources; from pathlib import Path; import sasmol; assert Path(sys.prefix).resolve() in Path(sasmol.__file__).resolve().parents; assert json.loads(resources.read_text('sasmol','_build_info.json')) == json.loads(Path(sys.argv[1]).read_text())"
    subprocess.run([sys.executable, "-I", "-c", script, str(expected)], check=True)
build_support.verify_installed_zazmol = verify_fixture
setup(name="sasmol", version="2.0", packages=["sasmol"],
      cmdclass=dict(build_support.build_commands(os.getcwd(), "zazmol", "2.0",
          "sasmol", build, build_py, sdist), install=build_support.verified_install_command(install)))
''')
    manifest = tmp_path / "served/release-manifest.json"
    manifest.parent.mkdir()
    manifest.write_bytes((app / "ui2/add/release-manifest.json").read_bytes())
    config = prefix / "etc/sassie/version-stack.json"
    config.parent.mkdir(parents=True)
    config.write_text(json.dumps(dict(schema_version=1, python_prefix=str(prefix),
                                     publisher=str(publisher), manifest=str(manifest))))
    result = subprocess.run([str(python), "setup.py", "install"], cwd=fixture, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    component = next(c for c in json.loads(manifest.read_text())["components"] if c["component_id"] == "zazmol")
    assert component["installed"]["version"] == "2.0"
    assert component["installed"]["status"] == "not_recorded"  # Fixture source has no Git receipt.
    before = manifest.read_bytes()
    result = subprocess.run([str(python), "setup.py", "install", "--root", str(tmp_path / "staging")],
                            cwd=fixture, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert manifest.read_bytes() == before
    assert "Staged/alternate installation" in result.stdout
    config.write_text(config.read_text().replace(str(manifest), str(tmp_path / "missing/release-manifest.json")))
    result = subprocess.run([str(python), "setup.py", "install"], cwd=fixture, text=True, capture_output=True)
    assert result.returncode != 0
    assert "Installed successfully; version-stack update failed" in result.stdout + result.stderr
    assert manifest.read_bytes() == before


def test_native_verifier_rejects_a_stale_stamp_and_python_extension_substitutes(tmp_path, monkeypatch, capfd):
    import json
    import venv
    import pytest
    support = load_build_support()
    prefix = tmp_path / "metadata environment"
    venv.EnvBuilder(with_pip=False).create(prefix)
    python = prefix / "bin/python"
    site = pathlib.Path(subprocess.check_output([str(python), "-I", "-c",
        "import sysconfig; print(sysconfig.get_path('purelib'))"], text=True).strip())
    package = site / "sasmol"
    package.mkdir()
    for name in ("__init__", "native_config", "_dcdio"):
        (package / (name + ".py")).write_text("")
    stamp = dict(schema_version=1, component_id="zazmol", version="2.0",
                 revision="a" * 40, modified=False, built_at="2026-10-15T12:00:00Z")
    (package / "_build_info.json").write_text(json.dumps(stamp))
    distribution = site / "sasmol-2.0.dist-info"
    distribution.mkdir()
    (distribution / "METADATA").write_text("Name: sasmol\nVersion: 2.0\n")
    expected = tmp_path / "verified-build.json"
    monkeypatch.setattr(support.sys, "executable", str(python))
    expected.write_text(json.dumps(dict(stamp, revision="b" * 40)))
    with pytest.raises(subprocess.CalledProcessError): support.verify_installed_zazmol(expected)
    assert "build stamp does not match" in capfd.readouterr().err
    expected.write_text(json.dumps(stamp))
    with pytest.raises(subprocess.CalledProcessError): support.verify_installed_zazmol(expected)
    assert "extension is not installed in this package" in capfd.readouterr().err
