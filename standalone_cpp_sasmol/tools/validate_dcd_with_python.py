#!/usr/bin/env python3
"""Validate a DCD with Python zazmol and print a compact summary."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dcd", type=Path)
    args = parser.parse_args()

    import sasmol.config as config
    import sasmol.system as system

    molecule = system.Molecule(0)
    molecule.read_dcd(str(args.dcd))
    coor = molecule.coor()
    total = numpy.sum(coor, dtype=config.CALC_DTYPE)

    print(f"dcd={args.dcd}")
    print(f"frames={coor.shape[0]}")
    print(f"atoms={coor.shape[1]}")
    print(f"coordinate_sum={total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
