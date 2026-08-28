#!/usr/bin/env python3
"""Fast convention checks for the rules in PROJECT_RULES.md.

Every check here corresponds to a real defect that reached the repository.
None of them build, mesh, or run the binary, so the whole file runs in
seconds and is meant to be run on every change:

    python3 -m test_system.test_conventions

Rule ids match PROJECT_RULES.md.
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

FAILURES: list[str] = []
PASSED = 0


def check(rule: str, name: str, ok: bool, detail: str = "") -> None:
    global PASSED
    if ok:
        PASSED += 1
        print(f"  ok    {rule}  {name}")
    else:
        FAILURES.append(f"{rule}  {name}: {detail}")
        print(f"  FAIL  {rule}  {name}: {detail}")


def skip(rule: str, name: str, why: str) -> None:
    print(f"  skip  {rule}  {name} ({why})")


# --- R1/R2: mesh index conventions ------------------------------------------

def find_mesh_cases() -> list[str]:
    out = []
    for d in ("work", "compset"):
        base = os.path.join(ROOT, d)
        if not os.path.isdir(base):
            continue
        for case in sorted(os.listdir(base)):
            for sub in (os.path.join(base, case), os.path.join(base, case, "fem_mesh_output")):
                if all(os.path.exists(os.path.join(sub, f))
                       for f in ("vert.txt", "fac.txt", "meshGeneralInfo.txt")):
                    out.append(sub)
                    break
    return out


def test_mesh_indexing() -> None:
    cases = find_mesh_cases()
    if not cases:
        skip("R1", "mesh cell area sane", "no meshed case on disk")
        return
    for c in cases[:3]:
        v = np.loadtxt(os.path.join(c, "vert.txt"))[:, :2]
        f = np.loadtxt(os.path.join(c, "fac.txt"), dtype=int)
        if f.min() == 1:
            f = f - 1
        p = v[f]
        x, y = p[:, :, 0], p[:, :, 1]
        area = 0.5 * np.abs(np.sum(x * np.roll(y, -1, 1) - np.roll(x, -1, 1) * y, axis=1)).sum()
        bbox = (v[:, 0].max() - v[:, 0].min()) * (v[:, 1].max() - v[:, 1].min())
        label = os.path.relpath(c, ROOT)
        # Scrambled indexing inflated this ratio to ~140x on a real case.
        check("R1", f"cell area <= bbox ({label})", area <= 1.05 * bbox,
              f"total cell area {area:.0f} vs bbox {bbox:.0f}")


def test_utilities_guard_index_base() -> None:
    for name in ("plotMeshNearFault.py", "checkMeshQuality.py", "plotMeshWithFaults.py"):
        p = os.path.join(SCRIPTS, name)
        if not os.path.exists(p):
            skip("R1", f"{name} guards index base", "not present")
            continue
        src = open(p).read()
        check("R1", f"{name} guards index base",
              "min() == 1" in src and "- 1" in src,
              "no `if min() == 1` guard before shifting indices")


def test_nsmp_not_filtered() -> None:
    p = os.path.join(SCRIPTS, "plotMeshWithFaults.py")
    if not os.path.exists(p):
        skip("R2", "nsmp sliced by nfn", "not present")
        return
    src = open(p).read()
    check("R2", "nsmp not filtered on id > 0", "ids > 0" not in src,
          "filtering ids > 0 drops node 0 under 0-based indexing")


# --- R3/R4: fault geometry --------------------------------------------------

def test_geometry_checker() -> None:
    try:
        from check_fault_geometry import check_geometry
    except Exception as e:  # noqa: BLE001
        check("R4", "check_fault_geometry importable", False, str(e))
        return

    dup = {"a": np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])}
    sev = [s for s, _ in check_geometry(dup, 400.0)]
    check("R3", "duplicate x is BLOCK", "BLOCK" in sev, f"got {sev}")

    shared = {"a": np.array([[0.0, 0.0], [10.0, 0.0]]),
              "b": np.array([[10.0, 0.0], [20.0, 0.0]])}
    msgs = " ".join(m for _, m in check_geometry(shared, 400.0))
    check("R4", "shared node detected", "SHARED NODE" in msgs, msgs or "no finding")

    short = {"a": np.array([[0.0, 0.0], [5.0, 0.0]])}
    msgs = " ".join(m for _, m in check_geometry(short, 400.0))
    check("R4", "short fault detected", "RESOLUTION" in msgs, msgs or "no finding")

    fine = {"a": np.column_stack([np.arange(0.0, 60.0, 0.5), np.zeros(120)]),
            "b": np.column_stack([np.arange(0.0, 60.0, 0.5), np.full(120, 30.0)])}
    check("R4", "clean geometry gives no findings", check_geometry(fine, 400.0) == [],
          str(check_geometry(fine, 400.0)))


# --- R19: compset parameters must be assigned to par ------------------------

def test_compset_params_use_par_prefix() -> None:
    """A bare `name = value` in user_defined_params.py is silently discarded.

    case.setup reads attributes off `par`, so an assignment without the
    prefix sets a module-level local that nothing ever reads. It fails
    silently and only bites when the value differs from the default.
    """
    import glob
    bad = []
    for path in sorted(glob.glob(os.path.join(ROOT, "compset", "*", "user_defined_params.py"))):
        for i, line in enumerate(open(path), start=1):
            stripped = line.strip()
            if (not stripped or stripped.startswith("#")
                    or stripped.startswith(("from ", "import ", "par.", "par ="))):
                continue
            if re.match(r"^[A-Za-z_][A-Za-z0-9_]*\s*=", stripped):
                bad.append(f"{os.path.relpath(path, ROOT)}:{i}: {stripped}")
    check("R19", "compset params all assigned to par", not bad,
          "; ".join(bad))


# --- R20: mesh export must not drop elements --------------------------------

def test_mesh_export_keeps_all_elements() -> None:
    """fac.txt must hold every 2D element the .msh does.

    meshGenLib writes only quads. Where gmsh cannot recombine a region it
    leaves triangles, and those are silently discarded, leaving holes in the
    FE mesh. Nodes bordering a hole lose part of their cell fan and their
    split node orphans. See issue #1.
    """
    import glob
    checked = 0
    for msh in sorted(glob.glob(os.path.join(ROOT, "work", "*", "fem_mesh_output",
                                             "eqdynaMesh.msh")))[:3]:
        fac = os.path.join(os.path.dirname(msh), "fac.txt")
        if not os.path.exists(fac):
            continue
        try:
            import meshio
            m = meshio.read(msh)
        except Exception as e:  # noqa: BLE001
            skip("R20", "mesh export keeps all elements", f"meshio: {e}")
            return
        n2d = sum(len(cb.data) for cb in m.cells if cb.type in ("quad", "triangle"))
        nfac = sum(1 for _ in open(fac))
        label = os.path.relpath(msh, ROOT)
        check("R20", f"no elements dropped ({label})", nfac == n2d,
              f"fac.txt has {nfac} rows, .msh has {n2d} 2D elements "
              f"({n2d - nfac} dropped)")
        checked += 1
    if not checked:
        skip("R20", "mesh export keeps all elements", "no meshed case with a .msh")


# --- R18: fault side classification -----------------------------------------

def test_side_classification_uses_normal() -> None:
    try:
        from meshGenLib import judgeElemDirect
    except Exception as e:  # noqa: BLE001
        check("R18", "judgeElemDirect importable", False, str(e))
        return
    # A steeply running fault: tangent mostly +y. Cells on either side both
    # have dy > 0, so the old global-dy test called both "above" and the
    # master node orphaned. With the tangent, they must land on opposite
    # sides.
    tang = (0.2, 1.0)
    left = judgeElemDirect((0.0, 0.0), (-0.5, 0.6), ftTangent=tang)
    right = judgeElemDirect((0.0, 0.0), (0.5, 0.6), ftTangent=tang)
    above = {1, 2}
    check("R18", "steep fault: opposite sides classified differently",
          (left in above) != (right in above),
          f"left->{left}, right->{right}; both on the same side orphans the master")
    # Gentle fault along +x must agree with the historical convention.
    flat = (1.0, 0.0)
    check("R18", "gentle fault agrees with global-dy convention",
          judgeElemDirect((0.0, 0.0), (0.1, 0.5), ftTangent=flat) in above
          and judgeElemDirect((0.0, 0.0), (0.1, -0.5), ftTangent=flat) not in above,
          "sign convention flipped for a gently sloping fault")


# --- R11: dropped-E exponent rule -------------------------------------------

def test_dropped_e_truncates() -> None:
    p = os.path.join(SCRIPTS, "paleo_site_stats.py")
    if not os.path.exists(p):
        skip("R11", "dropped-E truncates", "paleo_site_stats.py absent")
        return
    import importlib.util
    spec = importlib.util.spec_from_file_location("pss", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    conv = mod._RepairingConverter()
    got = conv(b"0.8384675-101")
    # MATLAB load() reads this as 0.8384675, NOT 0.8384675e-101 (verified
    # against R2024b). Reproducing published numbers requires truncation.
    check("R11", "dropped-E truncates like MATLAB", abs(got - 0.8384675) < 1e-12,
          f"got {got!r}, expected 0.8384675")
    check("R11", "well-formed exponent untouched",
          abs(conv(b"0.5000000E-101") - 5e-102) < 1e-115, "E-form mangled")


# --- R13: derived constants -------------------------------------------------

def test_observed_site_constants() -> None:
    try:
        from saf_result_utils import OBSERVED_EQDYNA_X_KM as obs
    except Exception as e:  # noqa: BLE001
        skip("R13", "observed-rate site x", f"import failed: {e}")
        return
    expect = np.array([-164.315662175388, 21.621033567916, 129.481611408865])
    check("R13", "observed-rate site x matches MATLAB chain",
          np.allclose(np.asarray(obs, float), expect, atol=1e-6),
          f"got {np.asarray(obs, float)}, expected {expect}")


# --- R12: magnitude convention ----------------------------------------------

def test_magnitude_constants_documented() -> None:
    p = os.path.join(SCRIPTS, "plotRuptureDynamics")
    src = open(p).read() if os.path.exists(p) else ""
    check("R12", "plotRuptureDynamics uses model rigidity",
          "rou * vs * vs" in src or "rou*vs*vs" in src, "rigidity expression changed")
    claude = os.path.join(ROOT, "CLAUDE.md")
    doc = open(claude).read() if os.path.exists(claude) else ""
    check("R12", "magnitude offset documented",
          "3500" in doc and "0.04" in doc,
          "CLAUDE.md must state the 0.04 Mw offset vs the paper's rigidity")


# --- R15/R16/R17: case I/O and run.sh ---------------------------------------

def test_scripts_search_araw() -> None:
    for name in ("plotRuptureDynamics", "paleo_site_stats.py", "saf_result_utils.py"):
        p = os.path.join(SCRIPTS, name)
        if not os.path.exists(p):
            skip("R15", f"{name} searches aRawSimuData", "not present")
            continue
        check("R15", f"{name} searches aRawSimuData",
              "aRawSimuData" in open(p).read(), "no aRawSimuData fallback")


def test_case_setup_run_sh() -> None:
    p = os.path.join(SCRIPTS, "case.setup")
    if not os.path.exists(p):
        skip("R16", "run.sh template", "case.setup absent")
        return
    src = open(p).read()
    # Strip trailing comments before matching: the C_mesh=2 line legitimately
    # carries "# keep binaryop: it is the restart state".
    rm_lines = [l.split("#")[0] for l in src.splitlines() if "rm -f totalop.txt" in l]
    check("R16", "run.sh does not delete binaryop",
          all("binaryop" not in l for l in rm_lines),
          "rm -f line still removes binaryop, breaking restart")
    check("R17", "run wrapper is detached", "setsid" in src,
          "no setsid; a process-group kill orphans the binary")
    check("R17", "run.sh emits Job finished", "Job finished" in src, "marker missing")


def main() -> None:
    print("EQdyna.2Dcycle convention checks (PROJECT_RULES.md)\n")
    for fn in (test_mesh_indexing, test_utilities_guard_index_base, test_nsmp_not_filtered,
               test_side_classification_uses_normal, test_compset_params_use_par_prefix,
               test_mesh_export_keeps_all_elements,
               test_geometry_checker, test_dropped_e_truncates, test_observed_site_constants,
               test_magnitude_constants_documented, test_scripts_search_araw,
               test_case_setup_run_sh):
        fn()
    print(f"\n{PASSED} passed, {len(FAILURES)} failed")
    for f in FAILURES:
        print(f"  {f}")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()


def test_no_hardcoded_fault_counts_in_fortran() -> None:
    """Fortran must never enumerate faults by hand -- loop over ntotft.

    Three separate bugs came from hand-written per-fault branches:

      interstress.f90  excluded fault tips and searched for nucleation on
                       faults 1-3 only, so gulang (5 faults) ran with ft4/ft5
                       unable to host an earthquake at all.
      faulting.f90     resolved the nucleation point with four if/elseif
                       branches and NO else, so a nucleation on fault 5+ left
                       ift0/xcoor0/ycoor0 uninitialised. The forced-rupture
                       patch was then never applied: the nucleating node crept
                       a few microns, stayed above failure, and every later
                       interseismic period collapsed to the 1-year floor. No
                       error, no warning -- just no earthquakes.

    Both were silent on 3-fault cases and only surfaced on 5- and 7-fault
    compsets. The mechanical guard is to ban the syntactic pattern: any
    expression summing three or more literal nfnode(<digit>) terms.
    """
    import re, glob
    pattern = re.compile(r"nfnode\(\d\)(\s*\+\s*nfnode\(\d\)){2,}")
    offenders = []
    for f in sorted(glob.glob(os.path.join(ROOT, "src", "*.f90"))):
        with open(f) as fh:
            for n, line in enumerate(fh, 1):
                if line.lstrip().startswith("!"):
                    continue
                if pattern.search(line):
                    offenders.append(f"{os.path.relpath(f, ROOT)}:{n}: {line.strip()[:90]}")
    assert not offenders, (
        "hardcoded per-fault branching in Fortran (loop over ntotft instead):\n  "
        + "\n  ".join(offenders))
