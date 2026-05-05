#!/usr/bin/env python3
"""Emit normalized JSON from Python zazmol's CHARMM topology parser.

This is an oracle helper for the standalone C++ port. It intentionally reports
Python's parsed shape, including string-valued charges and nested DELE records,
without trying to coerce the data into the eventual C++ model.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any


def repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[2]


def normalize(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: normalize(value[key]) for key in sorted(value)}
    if isinstance(value, tuple):
        return [normalize(item) for item in value]
    if isinstance(value, list):
        return [normalize(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Normalize Python CharmmTopology parser output as JSON."
    )
    parser.add_argument("topology_fixture", help="Path to a CHARMM topology fixture")
    args = parser.parse_args()

    root = repo_root()
    sys.path.insert(0, str(root / "src" / "python"))

    import sasmol.charmm_topology as charmm_topology

    fixture = pathlib.Path(args.topology_fixture).resolve()
    topology = charmm_topology.CharmmTopology()
    errors = topology.read_charmm_topology(
        topology_file_path=str(fixture.parent),
        topology_file_name=fixture.name,
    )

    payload = {
        "fixture": fixture.name,
        "errors": errors,
        "topology_info": normalize(topology.topology_info),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
