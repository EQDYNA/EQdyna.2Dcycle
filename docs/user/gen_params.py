#!/usr/bin/env python3
"""Generate the parameter table in parameters.md from scripts/defaultParameters.py.

defaultParameters.py is the class case.setup reads every default from, so a
reference generated from it cannot describe a knob that does not exist or keep
a stale default. Each top-level assignment in `class parameters` becomes one
row: name, default, and the comment on its own line.

Usage:
    python3 docs/user/gen_params.py          # rewrite the table in parameters.md
    python3 docs/user/gen_params.py --check  # exit 1 if parameters.md is stale
"""

from __future__ import annotations

import ast
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SOURCE = os.path.join(ROOT, "scripts", "defaultParameters.py")
TARGET = os.path.join(os.path.dirname(os.path.abspath(__file__)), "parameters.md")
BEGIN = ("<!-- BEGIN PARAMETER REFERENCE (generated from scripts/defaultParameters.py "
         "by docs/user/gen_params.py; do not edit by hand) -->")
END = "<!-- END PARAMETER REFERENCE -->"


def rows() -> list[tuple[str, str, str]]:
    src = open(SOURCE).read()
    lines = src.splitlines()
    tree = ast.parse(src)
    cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "parameters")
    out = []
    for node in cls.body:
        if not isinstance(node, ast.Assign):
            continue
        line = lines[node.lineno - 1]
        code, _, comment = line.partition("#")
        target = node.targets[0]
        names = ([e.id for e in target.elts] if isinstance(target, ast.Tuple)
                 else [target.id])
        default = ast.get_source_segment(src, node.value).strip()
        if isinstance(node.value, ast.Tuple) and len(names) == len(node.value.elts):
            defaults = [ast.get_source_segment(src, e).strip() for e in node.value.elts]
        else:
            defaults = [default] * len(names)
        for n, d in zip(names, defaults):
            if n == "exe":
                d = "`run_eqdyna2d_<VERSION>`"
            else:
                d = f"`{d}`"
            out.append((f"`{n}`", d, comment.strip().rstrip(".")))
    return out


def table() -> str:
    lines = ["| parameter | default | meaning |", "|---|---|---|"]
    for n, d, c in rows():
        lines.append(f"| {n} | {d} | {c.replace('|', '/')} |")
    return "\n".join(lines)


def main() -> int:
    page = open(TARGET).read()
    i, j = page.index(BEGIN) + len(BEGIN), page.index(END)
    new = page[:i] + "\n" + table() + "\n" + page[j:]
    if "--check" in sys.argv:
        if new != page:
            print("parameters.md is stale: run python3 docs/user/gen_params.py")
            return 1
        print("parameters.md is up to date")
        return 0
    open(TARGET, "w").write(new)
    print(f"wrote {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
