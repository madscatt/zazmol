"""Build-time provenance. Installed readers never consult a checkout."""

import datetime
import json
from pathlib import Path
import re
import subprocess
import sys
from setuptools.command.bdist_egg import bdist_egg


def build_info(root, component_id, version):
    root = Path(root).resolve()
    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(root), *args], stderr=subprocess.DEVNULL,
            universal_newlines=True).strip()
    revision, modified = None, False
    # A parent checkout is not evidence for a copied source tree.
    if (root / ".git").exists():
        revision = git("rev-parse", "--verify", "HEAD")
        modified = bool(git("status", "--porcelain", "--untracked-files=no"))
    elif (root / ".build-info.json").is_file():
        saved = json.loads((root / ".build-info.json").read_text())
        if (saved.get("schema_version") != 1 or
                saved.get("component_id") != component_id or
                saved.get("version") != version or
                type(saved.get("modified")) is not bool):
            raise ValueError("invalid source archive build information")
        revision, modified = saved.get("revision"), saved["modified"]
    if revision is not None and not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("invalid source revision")
    return {"schema_version": 1, "component_id": component_id,
            "version": version, "revision": revision, "modified": modified,
            "built_at": datetime.datetime.now(datetime.timezone.utc).replace(
                microsecond=0).isoformat()}


def write_info(path, info):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")


def build_commands(root, component_id, version, package, build, build_py, sdist):
    """Compose with each owner's native build commands, without changing them."""
    class StampedBuild(build):
        def run(self):
            target = Path(self.build_lib) / package / "_build_info.json"
            if target.exists():
                target.unlink()
            before = build_info(root, component_id, version)
            super().run()
            after = build_info(root, component_id, version)
            if any(before[k] != after[k] for k in ("revision", "modified")):
                raise RuntimeError("source changed during build; build again")
            write_info(target, after)

    class StampedBuildPy(build_py):
        def get_outputs(self, include_bytecode=1):
            outputs = super().get_outputs(include_bytecode)
            stamp = str(Path(self.build_lib) / package / "_build_info.json")
            return outputs if stamp in outputs else outputs + [stamp]

    class StampedEgg(bdist_egg):
        def run(self):
            # Older setuptools goes straight to build_py for eggs. Complete
            # the owner's entire native build before staging the stamped egg.
            self.run_command("build")
            super().run()

    class StampedSdist(sdist):
        def make_release_tree(self, base_dir, files):
            super().make_release_tree(base_dir, files)
            self.copy_file(str(Path(root) / "build_support.py"),
                           str(Path(base_dir) / "build_support.py"))
            write_info(Path(base_dir) / ".build-info.json",
                       build_info(root, component_id, version))

    return {"build": StampedBuild, "build_py": StampedBuildPy,
            "sdist": StampedSdist, "bdist_egg": StampedEgg}


def complete_install(component_id, expected_stamp):
    """Notify an explicitly bound instance only after owner verification."""
    config_path = Path(sys.prefix) / "etc/sassie/version-stack.json"
    try:
        if not config_path.exists() and not config_path.is_symlink():
            print("Installed build recorded; this Python environment has no bound website.")
            return
        config = json.loads(config_path.read_text())
        if (not isinstance(config, dict) or config.get("schema_version") != 1 or
                not isinstance(config.get("python_prefix"), str) or
                not Path(config["python_prefix"]).is_absolute() or
                Path(config["python_prefix"]).resolve() != Path(sys.prefix).resolve() or
                not isinstance(config.get("publisher"), str) or
                not Path(config["publisher"]).is_absolute() or
                not Path(config["publisher"]).is_file()):
            raise ValueError("invalid version-stack binding for this Python environment")
        subprocess.run([sys.executable, "-I", config["publisher"], "complete-install",
                        "--component", component_id, "--expected-stamp",
                        str(Path(expected_stamp).resolve()), "--verified"], check=True)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise RuntimeError("Installed successfully; version-stack update failed: " + str(error)) from error


def verify_installed_zazmol(expected_stamp):
    """Check the new installed package in a fresh, isolated interpreter."""
    script = r'''
import importlib
from importlib import machinery, metadata, resources
import json
from pathlib import Path
import sys
import sasmol
from sasmol import native_config

def within(path, parent):
    return path == parent or parent in path.parents

prefix = Path(sys.prefix).resolve()
package = Path(sasmol.__file__).resolve().parent
if not within(package, prefix):
    raise SystemExit("installed sasmol is outside this Python environment")
stamp = json.loads(resources.read_text("sasmol", "_build_info.json"))
if (stamp != json.loads(Path(sys.argv[1]).read_text()) or
        stamp.get("component_id") != "zazmol" or
        stamp.get("version") != metadata.version("sasmol")):
    raise SystemExit("installed ZazMol build stamp does not match this build/distribution")
for name in ("_dcdio", "_mask", "_view_vmd", "overlap", "matrix_math"):
    module = importlib.import_module("sasmol." + name)
    path = Path(module.__file__).resolve()
    if not within(path, package) or not any(str(path).endswith(s) for s in machinery.EXTENSION_SUFFIXES):
        raise SystemExit("ZazMol extension is not installed in this package: " + name)
    print("Installed native extension:", path)
if native_config.native_api_version() < 1 or not native_config.has_modern_dcd_reader_api():
    raise SystemExit("installed ZazMol lacks the required modern native API")
paths = ([Path(p) / "sasmol/molecule.hpp" for p in native_config.include_dirs()] +
         [Path(p) for p in native_config.extra_objects()])
for path in paths:
    if not within(path.resolve(), package) or not path.is_file() or not path.stat().st_size:
        raise SystemExit("installed ZazMol native file is missing or outside the package: " + str(path))
    print("Installed native file:", path)
for directory in native_config.library_dirs():
    if not within(Path(directory).resolve(), package) or not Path(directory).is_dir():
        raise SystemExit("installed ZazMol native library directory is invalid")
print("Installed ZazMol verified:", stamp)
'''
    subprocess.run([sys.executable, "-I", "-c", script, str(expected_stamp)], check=True)


def verified_install_command(install):
    """Keep packaging/staging separate from a normal environment install."""
    class VerifiedInstall(install):
        def run(self):
            import sysconfig
            environment_install = (
                not (self.root or self.home or self.user or self.dry_run) and
                self.prefix is not None and
                Path(self.prefix).resolve() == Path(sys.prefix).resolve() and
                Path(self.install_lib).resolve() in {
                    Path(sysconfig.get_path(key)).resolve() for key in ("purelib", "platlib")})
            super().run()
            if not environment_install:
                print("Staged/alternate installation; instance version stack unchanged.")
                return
            build = self.get_finalized_command("build")
            expected_stamp = Path(build.build_lib) / "sasmol/_build_info.json"
            verify_installed_zazmol(expected_stamp)
            complete_install("zazmol", expected_stamp)
    return VerifiedInstall
