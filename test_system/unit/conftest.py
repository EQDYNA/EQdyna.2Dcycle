"""Unit tier: fast, pure-Python tests of the tools in scripts/ and docs/user/.

No Fortran build, no meshing, no simulation; the whole tier runs in seconds.
"""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
for p in (os.path.join(ROOT, "scripts"), os.path.join(ROOT, "docs", "user")):
    if p not in sys.path:
        sys.path.insert(0, p)
