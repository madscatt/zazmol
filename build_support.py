"""Build-time provenance. Installed readers never consult a checkout."""

import datetime
import json
from pathlib import Path
import re
import subprocess
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
