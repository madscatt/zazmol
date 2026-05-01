#!/usr/bin/env python3
"""Emit a read-only Markdown survey of Python zazmol public methods."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path


def public_name(name: str) -> bool:
    return not name.startswith("_")


def signature_text(node: ast.FunctionDef) -> str:
    pieces: list[str] = []
    defaults = [None] * (len(node.args.args) - len(node.args.defaults))
    defaults.extend(node.args.defaults)

    for arg, default in zip(node.args.args, defaults):
        if default is None:
            pieces.append(arg.arg)
        else:
            pieces.append(f"{arg.arg}=...")

    if node.args.vararg:
        pieces.append(f"*{node.args.vararg.arg}")

    for arg in node.args.kwonlyargs:
        pieces.append(f"{arg.arg}=...")

    if node.args.kwarg:
        pieces.append(f"**{node.args.kwarg.arg}")

    return f"{node.name}({', '.join(pieces)})"


def module_entries(path: Path) -> list[tuple[str, str, int]]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    entries: list[tuple[str, str, int]] = []

    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and public_name(node.name):
            entries.append(("function", signature_text(node), node.lineno))
        elif isinstance(node, ast.ClassDef) and public_name(node.name):
            entries.append(("class", node.name, node.lineno))
            for child in node.body:
                if isinstance(child, ast.FunctionDef) and public_name(child.name):
                    entries.append(
                        ("method", f"{node.name}.{signature_text(child)}", child.lineno)
                    )

    return entries


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--python-src",
        type=Path,
        default=Path("src/python"),
        help="Path to the Python zazmol source directory",
    )
    args = parser.parse_args()

    python_src = args.python_src.resolve()
    if not python_src.exists():
        raise SystemExit(f"Python source directory not found: {python_src}")

    print("# Python zazmol Surface Survey")
    print()
    print(f"Source: `{python_src}`")
    print()
    print("| Module | Kind | Name | Line | C++ Status | Notes |")
    print("| --- | --- | --- | ---: | --- | --- |")

    for path in sorted(python_src.glob("*.py")):
        if path.name.startswith("__"):
            continue
        module = path.stem
        for kind, name, line in module_entries(path):
            print(f"| `{module}` | {kind} | `{name}` | {line} | question | |")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
