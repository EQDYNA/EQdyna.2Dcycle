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
    # The generated run.sh used to re-run meshgen.py on every launch and then
    # unconditionally `cp fem_mesh_output/<file> .`. That silently discarded
    # per-node loading patched on top of the mesh (apply_strain_loading.py
    # rewrites nsmpGeoPhys.txt in place), and the run proceeded with default
    # uniform loading with nothing in the log to say so.
    check("R18", "run.sh skips meshing when a mesh exists",
          "FORCE_MESH" in src,
          "meshgen.py runs unconditionally; a re-mesh discards per-node loading")
    check("R18", "run.sh does not clobber newer mesh working copies",
          '-nt' in src and 'not overwriting' in src,
          "mesh files are copied unconditionally over locally patched versions")


# --- R24/R26/R27: release process -------------------------------------------

def _read_version() -> str:
    p = os.path.join(ROOT, "VERSION")
    with open(p) as fh:
        return fh.read().strip()


def test_version_file_is_single_semver_line() -> None:
    """R27: VERSION is the single source of the release version.

    A single line of the form MAJOR.MINOR.PATCH, optionally with a
    '-<prerelease>' suffix (e.g. 2.0.7-rc8). Anything else (blank, multiple
    lines, trailing garbage) means src/makefile's `$(shell cat VERSION)`
    produces a binary name nobody asked for.
    """
    p = os.path.join(ROOT, "VERSION")
    if not os.path.exists(p):
        check("R27", "VERSION file exists", False, "root VERSION file missing")
        return
    raw = open(p).read()
    lines = raw.splitlines()
    check("R27", "VERSION is exactly one line", len(lines) == 1,
          f"got {len(lines)} lines: {lines!r}")
    if lines:
        check("R27", "VERSION is a semver (+ optional -prerelease)",
              re.fullmatch(r"\d+\.\d+\.\d+(-[0-9A-Za-z.]+)?", lines[0]) is not None,
              f"got {lines[0]!r}")


def extract_changelog_section(version: str, changelog_path: str) -> str:
    """The correct CHANGELOG section extraction (R26).

    The naive awk range /^## \\[$VERSION\\]/,/^## \\[/ closes on the same
    line it opens, because the start heading also matches the end pattern.
    Every tag through v2.0.7-rc7 shipped with an empty message because of
    this. Skip the heading line itself, then stop at the *next* heading.
    """
    hdr = f"## [{version}]"
    lines_out = []
    found = False
    for line in open(changelog_path):
        if not found:
            if line.startswith(hdr):
                found = True
            continue
        if line.startswith("## ["):
            break
        lines_out.append(line)
    return "".join(lines_out)


def test_changelog_has_body_for_version() -> None:
    """R24: CHANGELOG must have a real body under '## [<VERSION>]'."""
    changelog = os.path.join(ROOT, "CHANGELOG.md")
    if not os.path.exists(changelog):
        check("R24", "CHANGELOG.md has a section for VERSION", False, "CHANGELOG.md missing")
        return
    version = _read_version()
    src = open(changelog).read()
    hdr_present = re.search(rf"^## \[{re.escape(version)}\]", src, re.MULTILINE) is not None
    check("R24", f"CHANGELOG.md has a '## [{version}]' heading", hdr_present,
          "no heading for the current VERSION -- rename '## [Unreleased]' and add today's date")
    if not hdr_present:
        return
    body = extract_changelog_section(version, changelog)
    check("R24", f"CHANGELOG.md '## [{version}]' has a non-empty body",
          body.strip() != "",
          "heading exists but body is blank -- an empty release note")


def test_changelog_section_extraction_returns_text() -> None:
    """R26: the section-extraction a human or script uses must not be the
    self-closing awk range that produced an empty tag message on every
    release through v2.0.7-rc7."""
    changelog = os.path.join(ROOT, "CHANGELOG.md")
    if not os.path.exists(changelog):
        skip("R26", "changelog extraction returns text", "CHANGELOG.md missing")
        return
    version = _read_version()
    body = extract_changelog_section(version, changelog)
    check("R26", "extract_changelog_section(VERSION) returns non-empty text",
          body.strip() != "",
          f"extraction for {version!r} returned nothing -- this is exactly the "
          "self-closing-range bug if it recurs")


def test_run_sh_sets_omp_threads() -> None:
    """R30: results reproduce only at a fixed thread count, so run.sh must pin one.

    qdct3 and hrglss accumulate into brhs under $OMP ATOMIC, so the summation
    order depends on thread scheduling and the result differs in the last bit
    between thread counts. Earthquake sequences are chaotic, so that amplifies
    to a different catalogue within about four cycles. Verified on v2.1.0:
    at OMP_NUM_THREADS=1 both paper.saf.A (C_mesh=2) and saf.gmsh.lite
    (C_mesh=3) reproduce their stored totalop.txt1 exactly (max abs diff 0.0);
    at 2 threads the same binary diverges by cycle 5.

    If run.sh inherited the caller's OMP_NUM_THREADS, the thread count -- and
    so reproducibility -- would depend on unrecorded shell state.
    """
    p = os.path.join(SCRIPTS, "case.setup")
    if not os.path.exists(p):
        skip("R30", "run.sh pins thread count", "case.setup absent")
        return
    src = open(p).read()
    check("R30", "run.sh sets OMP_NUM_THREADS explicitly",
          "OMP_NUM_THREADS" in src,
          "generated run.sh inherits the caller's thread count; runs become "
          "irreproducible against stored references")


def test_every_check_is_registered() -> None:
    """Every test_* function must appear in main()'s run list.

    Twice now a check has been added to this file and left out of that tuple,
    so `python3 test_system/test_conventions.py` silently never ran it and only
    pytest's auto-discovery did -- a guard that looks installed but is not.
    It happened to test_no_hardcoded_fault_counts_in_fortran (the guard against
    the faulting.f90 / interstress.f90 fault-count bugs) and again to
    test_run_sh_sets_omp_threads.
    """
    import inspect
    src = inspect.getsource(main)
    here = sys.modules[__name__]
    missing = [n for n in dir(here)
               if n.startswith("test_") and n != "test_every_check_is_registered"
               and callable(getattr(here, n)) and n not in src]
    check("R21", "every test_* function is registered in main()",
          not missing,
          "unregistered, so the plain runner skips them: " + ", ".join(missing))


def test_fig4_refuses_empty_output() -> None:
    """R29: plot_event_slips_overtime_fig4.py must not ship an empty figure.

    It saved axes, fault traces and a scale bar with no events on them and
    exited 0 -- once when SAF-tuned thresholds filtered out every event of a
    sub-metre catalogue, once when the requested window started past the last
    event. Both looked like finished plots.
    """
    p = os.path.join(SCRIPTS, "plot_event_slips_overtime_fig4.py")
    if not os.path.exists(p):
        skip("R29", "fig4 refuses empty output", "script absent")
        return
    src = open(p).read()
    check("R29", "fig4 counts the events it draws",
          "shown == 0" in src or "if empty" in src,
          "no empty-output detection")
    check("R29", "fig4 exits non-zero when nothing is drawn",
          "sys.exit(1)" in src,
          "an empty figure still exits 0 and reads as a finished result")


def main() -> None:
    print("EQdyna.2Dcycle convention checks (PROJECT_RULES.md)\n")
    for fn in (test_mesh_indexing, test_utilities_guard_index_base, test_nsmp_not_filtered,
               test_side_classification_uses_normal, test_compset_params_use_par_prefix,
               test_mesh_export_keeps_all_elements,
               test_geometry_checker, test_dropped_e_truncates, test_observed_site_constants,
               test_magnitude_constants_documented, test_scripts_search_araw,
               test_case_setup_run_sh, test_no_hardcoded_fault_counts_in_fortran,
               test_version_file_is_single_semver_line, test_changelog_has_body_for_version,
               test_changelog_section_extraction_returns_text,
               test_run_sh_sets_omp_threads, test_fig4_refuses_empty_output,
               test_every_check_is_registered):
        fn()
    print(f"\n{PASSED} passed, {len(FAILURES)} failed")
    for f in FAILURES:
        print(f"  {f}")
    sys.exit(1 if FAILURES else 0)


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
    import glob
    pattern = re.compile(r"nfnode\(\d\)(\s*\+\s*nfnode\(\d\)){2,}")
    offenders = []
    for f in sorted(glob.glob(os.path.join(ROOT, "src", "*.f90"))):
        with open(f) as fh:
            for n, line in enumerate(fh, 1):
                if line.lstrip().startswith("!"):
                    continue
                if pattern.search(line):
                    offenders.append(f"{os.path.relpath(f, ROOT)}:{n}: {line.strip()[:90]}")
    check("R28", "no hardcoded per-fault branching in Fortran", not offenders,
          "loop over ntotft instead:\n  " + "\n  ".join(offenders))


if __name__ == "__main__":
    main()
