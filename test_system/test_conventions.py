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

HERE = os.path.dirname(os.path.abspath(__file__))
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


def test_stf_format_versions_match() -> None:
    """R31: stf.bin writer and reader must agree on the format version."""
    import re
    w = os.path.join(ROOT, "src", "stf_output.f90")
    r = os.path.join(SCRIPTS, "stf_read.py")
    if not (os.path.exists(w) and os.path.exists(r)):
        skip("R31", "stf writer/reader versions match", "files absent")
        return
    ws = open(w).read()
    rs = open(r).read()
    mw = re.search(r"write\(u\) int\((\d+),4\), int\(icstart", ws)
    mr = re.search(r"if self\.version != (\d+)", rs)
    check("R31", "stf.bin writer and reader carry the same format version",
          bool(mw and mr and mw.group(1) == mr.group(1)),
          f"writer {mw.group(1) if mw else '?'} vs reader {mr.group(1) if mr else '?'}")
    check("R31", "stf.bin magic matches", "'EQDYNSTF'" in ws and 'b"EQDYNSTF"' in rs,
          "magic string differs between writer and reader")


def test_user_guide() -> None:
    """The user guide must build from current code and stay user-facing.

    docs/user/parameters.md is generated from scripts/defaultParameters.py; a
    stale table documents defaults the code no longer has. User pages must not
    cite internal material (rule ids, the rule book, agent names) that a reader
    outside the project cannot follow.
    """
    import re, subprocess, glob
    d = os.path.join(ROOT, "docs", "user")
    if not os.path.isdir(d):
        skip("R32", "user guide", "docs/user absent")
        return
    r = subprocess.run([sys.executable, os.path.join(d, "gen_params.py"), "--check"],
                       capture_output=True, text=True)
    check("R32", "parameters.md matches defaultParameters.py", r.returncode == 0,
          "stale: run python3 docs/user/gen_params.py")
    bad = re.compile(r"\bR\d{1,2}\b|PROJECT_RULES|\b(?:zofia|victor|mira|kai|haruto|wei-lin|"
                     r"lars|iris|sophia|nadia)\b|board row", re.I)
    hits = []
    for f in sorted(glob.glob(os.path.join(d, "*.md"))):
        for n, line in enumerate(open(f), 1):
            if "BEGIN PARAMETER REFERENCE" in line:
                continue
            if bad.search(line):
                hits.append(f"{os.path.basename(f)}:{n}")
    check("R32", "user guide has no internal references", not hits, ", ".join(hits[:5]))


def test_no_tracked_bytecode() -> None:
    """R33: compiled Python never enters the repository.

    The ignore rules once named __pycache__ per directory; a new test directory
    was not covered and six .pyc files were committed with it.
    """
    import subprocess
    r = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        skip("R33", "no tracked bytecode", "not a git checkout")
        return
    bad = [f for f in r.stdout.splitlines() if f.endswith((".pyc", ".pyo")) or "__pycache__/" in f]
    check("R33", "no .pyc / __pycache__ tracked", not bad, ", ".join(bad[:5]))


def test_root_whitelist() -> None:
    """R35: the repo root is a whitelist, not a preference."""
    import subprocess
    r = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        skip("R35", "root whitelist", "not a git checkout")
        return
    allowed_files = {
        ".gitignore", "CHANGELOG.md", "CLAUDE.md", "LICENSE", "PROJECT_RULES.md",
        "README.md", "VERSION", "example_workflow.sh", "install.sh",
        "resume_claude.sh", "PATHWAY_FORWARD.md",
    }
    allowed_dirs = {".github", "compset", "docs", "misc", "scripts", "src", "test_system"}
    offenders = []
    for f in r.stdout.splitlines():
        if "/" not in f:
            if f not in allowed_files:
                offenders.append(f)
        else:
            top = f.split("/", 1)[0]
            if top not in allowed_dirs:
                offenders.append(f)
    check("R35", "tracked root matches the whitelist", not offenders,
          "unexpected root entries: " + ", ".join(sorted(set(
              o if "/" not in o else o.split('/', 1)[0] for o in offenders))[:10]))


def test_root_file_size() -> None:
    """R35a: no tracked file over 5 MB."""
    import subprocess
    r = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True)
    if r.returncode != 0:
        skip("R35a", "no tracked file over 5 MB", "not a git checkout")
        return
    files = [f for f in r.stdout.decode().split("\0") if f]
    limit = 5 * 1024 * 1024
    offenders = []
    for f in files:
        path = os.path.join(ROOT, f)
        try:
            size = os.path.getsize(path)
        except OSError:
            continue
        if size > limit:
            offenders.append(f"{f} ({size // (1024*1024)} MB)")
    check("R35a", "no tracked file over 5 MB", not offenders, ", ".join(offenders[:10]))


def test_repo_tidy() -> None:
    """R35a: report-only tidy pass -- stale worktrees, merged branches, uncited runs/.

    Never fails the suite; prints findings for the owner to act on.
    """
    import subprocess
    notes = []

    wt = subprocess.run(["git", "worktree", "list", "--porcelain"], cwd=ROOT,
                         capture_output=True, text=True)
    worktrees = [l.split(" ", 1)[1] for l in wt.stdout.splitlines() if l.startswith("worktree ")]
    stale = [w for w in worktrees if os.path.normpath(w) != os.path.normpath(ROOT)
             and not os.path.isdir(w)]
    if stale:
        notes.append(f"stale worktree entries: {', '.join(stale)}")

    # `git branch --merged` misses squash merges (every PR here is squashed).
    # A branch is merged when its whole diff, squashed onto its merge base,
    # already has an equivalent patch in main (`git cherry` prints "-").
    def git(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    refs = git("for-each-ref", "--format=%(refname:short)", "refs/heads", "refs/remotes/origin").split()
    merged = []
    for ref in refs:
        if ref in ("main", "origin/main", "origin/HEAD", "origin"):
            continue
        base = git("merge-base", "main", ref)
        if not base:
            continue
        if git("rev-parse", ref) == base:
            merged.append(ref)
            continue
        squash = git("commit-tree", f"{ref}^{{tree}}", "-p", base, "-m", "tidy-probe")
        if squash and git("cherry", "main", squash).startswith("-"):
            merged.append(ref)
    if merged:
        notes.append(f"branches already merged into main (incl. squash): {', '.join(merged)}")

    runs_dir = os.path.join(ROOT, "runs")
    if os.path.isdir(runs_dir):
        entries = sorted(os.listdir(runs_dir))
        if entries:
            notes.append(f"runs/ entries present (verify still cited): {', '.join(entries[:10])}")

    if notes:
        print("  note  R35a  repo tidy (report-only): " + "; ".join(notes))
    else:
        print("  note  R35a  repo tidy (report-only): nothing to report")


def test_no_unmarked_silent_fallback() -> None:
    """R36: `|| true` needs a stated reason on the same line.

    install.sh and make_results_bundle.sh both swallowed real cp/rm failures
    this way; fixed by dropping the unneeded `|| true` where `set -e` was
    already enough. Where it is genuine ("this glob may legitimately not
    match yet"), it must carry an inline `# intentional: <why>` marker.
    """
    import glob
    sh_files = sorted(glob.glob(os.path.join(ROOT, "*.sh")) +
                       glob.glob(os.path.join(ROOT, "scripts", "*.sh")))
    # The run.sh template is Python source (scripts/case.setup) that emits
    # shell text via f.write(...); check the emitted text, not the .py file.
    offenders = []
    for path in sh_files:
        lines = open(path).readlines()
        for i, line in enumerate(lines, 1):
            if "|| true" not in line:
                continue
            prev = lines[i - 2] if i >= 2 else ""
            if "intentional:" not in line and "intentional:" not in prev:
                offenders.append(f"{os.path.relpath(path, ROOT)}:{i}")
    case_setup = os.path.join(SCRIPTS, "case.setup")
    if os.path.exists(case_setup):
        for i, line in enumerate(open(case_setup), 1):
            if "|| true" in line and "intentional:" not in line:
                offenders.append(f"scripts/case.setup:{i} (emitted into run.sh)")
    check("R36", "no unmarked `|| true`", not offenders, ", ".join(offenders[:10]))


def test_status_board_exists() -> None:
    """R37: PATHWAY_FORWARD.md exists, parses, and every row matches the header width."""
    p = os.path.join(ROOT, "PATHWAY_FORWARD.md")
    if not os.path.exists(p):
        check("R37", "PATHWAY_FORWARD.md exists", False, "missing at repo root")
        return
    rows = [l for l in open(p) if l.strip().startswith("|")]
    tables = []
    cur = []
    for l in rows:
        cur.append(l)
    # Group contiguous pipe-rows into tables (separated by non-pipe lines already excluded)
    if rows:
        header_cols = rows[0].count("|")
        bad = [f"line with {l.count('|')} cols (header has {header_cols}): {l.strip()[:60]}"
               for l in rows if l.count("|") != header_cols and not set(l.strip()) <= set("|-: ")]
        check("R37", "every board row matches the header's column count", not bad,
              "; ".join(bad[:5]))
    else:
        check("R37", "PATHWAY_FORWARD.md has at least one table", False, "no '|' rows found")


def test_no_rival_status_board() -> None:
    """R37: one board, everywhere -- no TODO/STATUS/BACKLOG/ROADMAP file beside it."""
    import subprocess, re as _re
    r = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        skip("R37", "no rival status board", "not a git checkout")
        return
    pat = _re.compile(r"(TODO|STATUS|ROADMAP|BACKLOG|TASKS|PLAN)\.(md|txt)$", _re.IGNORECASE)
    bad = [f for f in r.stdout.splitlines() if pat.search(f)]
    check("R37", "no rival status board tracked", not bad, ", ".join(bad[:5]))


def test_feglobal_malformed_lines() -> None:
	"""R39a — Malformed STF/slipSense lines in FE_Global.txt hard-fail, not silent-default.
	
	PR #18 (00e7865) fixed an EOF escalation bug where the binary would abort 
	when FE_Global.txt lacked both STF and slipSense lines (old pre-2.3.0 format).
	The fix made it skip the slipSense read if STF read failed, but didn't 
	distinguish between "line missing" (ios < 0, EOF) and "line malformed" 
	(ios > 0, parse error). A malformed line should hard-fail loudly, not 
	silently default to slipSense=+1, since that would silently flip fault 
	sense without warning.
	
	This test verifies:
	- Malformed STF line: hard-fails with error message
	- Valid STF + malformed slipSense: hard-fails with error message  
	- Missing both STF and slipSense (old format): silently defaults to +1
	"""
	import subprocess
	import shutil
	import tempfile
	
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	version = open(os.path.join(root, "VERSION")).read().strip()
	exe = os.path.join(root, "bin", f"run_eqdyna2d_{version}")
	
	if not os.path.exists(exe):
		skip("R39a", "FE_Global malformed lines", "binary not built")
		return
	
	# Find a case to template from (any compset with FE_Global.txt already on
	# disk); failing that, generate one with case.setup -- a fresh checkout or
	# CI run never has a pre-existing FE_Global.txt (it's a generated file,
	# excluded from compset/ on purpose), so relying on one being present would
	# make this check skip, not run, on every clean environment.
	case_template = None
	for base in ("compset", "work"):
		base_path = os.path.join(root, base)
		if not os.path.isdir(base_path):
			continue
		for case in os.listdir(base_path):
			fg = os.path.join(base_path, case, "FE_Global.txt")
			if os.path.exists(fg):
				case_template = os.path.join(base_path, case)
				break
		if case_template:
			break

	generated_case = None
	if not case_template:
		import subprocess as _sp
		import tempfile as _tf
		generated_case = _tf.mkdtemp(prefix="r39a_case_")
		create_newcase = os.path.join(root, "scripts", "create.newcase")
		case_env = os.environ.copy()
		case_env["EQDYNA2DCYCLEROOT"] = root
		r = _sp.run(["python3", create_newcase, "--work_dir", generated_case,
		             "--compset", "gulang.gmsh.lite", "--force"],
		            cwd=root, capture_output=True, text=True, timeout=60, env=case_env)
		if r.returncode == 0:
			r2 = _sp.run(["python3", "case.setup"], cwd=generated_case,
			              capture_output=True, text=True, timeout=60, env=case_env)
			if r2.returncode == 0 and os.path.exists(os.path.join(generated_case, "FE_Global.txt")):
				case_template = generated_case
		if not case_template:
			shutil.rmtree(generated_case, ignore_errors=True)
			skip("R39a", "FE_Global malformed lines",
			     f"could not generate a case to template from: {r.stderr[:200] or (r2.stderr[:200] if 'r2' in dir() else '')}")
			return

	fails = []
	
	# Test 1: Malformed STF line (present but unparseable)
	case1 = tempfile.mkdtemp(prefix="test_malformed_stf_")
	try:
		shutil.copytree(case_template, case1, dirs_exist_ok=True)
		fg = os.path.join(case1, "FE_Global.txt")
		with open(fg, 'r') as f:
			lines = f.read().splitlines()
		# Replace STF line (usually line 27, 0-indexed) with garbage
		for i, line in enumerate(lines):
			if line and not line.isspace():
				# Find the STF line (the one that looks like "0 1 0.001" or similar)
				if i > 25 and i < 30:  # STF is near the end, after plotmesh
					try:
						parts = line.split()
						if len(parts) == 3:
							float(parts[0]); float(parts[1]); float(parts[2])
							lines[i] = "garbage not a number xyz"
							break
					except (ValueError, IndexError):
						pass
		with open(fg, 'w') as f:
			f.write("\n".join(lines) + "\n")
		
		shutil.copy2(exe, case1)
		env = os.environ.copy()
		env['OMP_NUM_THREADS'] = '1'
		r = subprocess.run([os.path.join(case1, os.path.basename(exe))], cwd=case1, 
		                   capture_output=True, text=True, timeout=30)
		
		if r.returncode != 0 and "error: could not parse STF line" in r.stdout:
			pass  # Good
		else:
			fails.append("malformed STF line should hard-fail")
	except Exception as e:
		fails.append(f"test 1 exception: {e}")
	finally:
		shutil.rmtree(case1, ignore_errors=True)
	
	# Test 2: Valid STF + Malformed slipSense line
	case2 = tempfile.mkdtemp(prefix="test_malformed_slip_")
	try:
		shutil.copytree(case_template, case2, dirs_exist_ok=True)
		fg = os.path.join(case2, "FE_Global.txt")
		with open(fg, 'r') as f:
			lines = f.read().splitlines()
		# Replace slipSense line (usually last line) with garbage
		lines[-1] = "not a valid slipsense"
		with open(fg, 'w') as f:
			f.write("\n".join(lines) + "\n")
		
		shutil.copy2(exe, case2)
		env = os.environ.copy()
		env['OMP_NUM_THREADS'] = '1'
		r = subprocess.run([os.path.join(case2, os.path.basename(exe))], cwd=case2, 
		                   capture_output=True, text=True, timeout=30)
		
		if r.returncode != 0 and "error: could not parse slipSense" in r.stdout:
			pass  # Good
		else:
			fails.append("malformed slipSense line should hard-fail")
	except Exception as e:
		fails.append(f"test 2 exception: {e}")
	finally:
		shutil.rmtree(case2, ignore_errors=True)
	
	# Test 3: Missing both STF and slipSense (old format)
	case3 = tempfile.mkdtemp(prefix="test_missing_stf_")
	try:
		shutil.copytree(case_template, case3, dirs_exist_ok=True)
		fg = os.path.join(case3, "FE_Global.txt")
		with open(fg, 'r') as f:
			lines = f.read().splitlines()
		# Truncate to remove STF and slipSense lines (pre-2.3.0 format)
		lines = lines[:27]  # Keep only lines up to plotmesh
		with open(fg, 'w') as f:
			f.write("\n".join(lines) + "\n")
		
		shutil.copy2(exe, case3)
		env = os.environ.copy()
		env['OMP_NUM_THREADS'] = '1'
		r = subprocess.run([os.path.join(case3, os.path.basename(exe))], cwd=case3, 
		                   capture_output=True, text=True, timeout=30)
		
		# Should silently default, not error
		if r.returncode == 0:
			pass  # Good - PR #18's original fix should still work
		else:
			fails.append("missing STF/slipSense should silently default (PR #18 regression)")
	except Exception as e:
		fails.append(f"test 3 exception: {e}")
	finally:
		shutil.rmtree(case3, ignore_errors=True)
	
	
	# Test 4: REGRESSION — STF line deleted, only lone -1 left
	# PR #21 regression: when STF line is genuinely missing but a slipSense-shaped
	# value (-1) is left as the sole trailing line, the OLD (broken) code would:
	#   1. Read the entire "-1" line looking for (stfOn, stfEvery, stfVmin)
	#   2. Get ios < 0 (EOF after finding -1 in stfOn position, missing stfEvery/stfVmin)
	#   3. Interpret as "STF line missing" and silently default
	#   4. Skip the slipSense read entirely
	#   5. Leave slipSense as its uninitialized value (+1 in Fortran)
	# This silently flipped left-lateral faults to right-lateral.
	# THE FIX: distinguish line presence (character read) from value parsing (internal read).
	# A line present but failing to parse must hard-fail, not default.
	case4 = tempfile.mkdtemp(prefix="test_stf_deleted_lone_")
	try:
		shutil.copytree(case_template, case4, dirs_exist_ok=True)
		fg = os.path.join(case4, "FE_Global.txt")
		with open(fg, 'r') as f:
			lines_fg = f.read().splitlines()
		# Remove the STF line (line 27, 0-indexed), leaving only the next line (-1)
		# This simulates a user who deleted the STF line but left a slipSense value
		# The file now has only 27 complete data lines + the lone "-1"
		lines_fg = lines_fg[:27]  # Remove lines 28 and 29
		lines_fg.append("-1")   # Add back just the "-1" as the new line 28
		with open(fg, 'w') as f:
			f.write("\n".join(lines_fg) + "\n")
		
		shutil.copy2(exe, case4)
		env = os.environ.copy()
		env['OMP_NUM_THREADS'] = '1'
		r = subprocess.run([os.path.join(case4, os.path.basename(exe))], cwd=case4, 
		                   capture_output=True, text=True, timeout=30)
		
		# The lone "-1" cannot parse as (stfOn, stfEvery, stfVmin) — it's missing
		# stfEvery and stfVmin. The internal read must fail (iosParse != 0), and
		# we must hard-fail, not silently default.
		if r.returncode != 0 and "error: could not parse STF line" in r.stdout:
			pass  # Good — rejected the truncated STF line
		else:
			fails.append("STF line deleted with lone -1: should hard-fail, not silent-default")
	except Exception as e:
		fails.append(f"test 4 exception: {e}")
	finally:
		shutil.rmtree(case4, ignore_errors=True)
	
	check("R39a", "FE_Global malformed lines hard-fail, missing lines silent-default",
	      not fails, "; ".join(fails[:3]))
	if generated_case:
		shutil.rmtree(generated_case, ignore_errors=True)


# --- R10: every bug gets a regression test ----------------------------------

def test_fe_global_pre230_no_crash() -> None:
    """A genuinely pre-2.3.0 FE_Global.txt (no STF or slipSense trailer at
    all) must not crash readglobal, and slipSense must silently default to
    +1 as documented.

    PR #17 (35a6549) made a malformed slipSense line a hard error. It did
    not account for gfortran escalating a second consecutive read past
    end-of-file from a negative (EOF) iostat to a positive (error) iostat:
    any FE_Global.txt missing *both* trailers -- i.e. every file written
    before this parameter existed -- hit the STF read's correct EOF
    default, then the slipSense read's new hard-error path, and the binary
    aborted with "could not parse slipSense line in FE_Global.txt". PR #18
    (00e7865) fixed it by skipping the slipSense read entirely once the STF
    read had already failed. CI never caught the original regression
    because every smoke/regression fixture has a current-format
    FE_Global.txt; only a hand-run case with an old file crashed (caught by
    e2e re-run before tagging v2.3.1, not by CI). This guards that gap: the
    fixture is a real FE_Global.txt from before v2.3.0 (git-archived from
    commit 2211914, "Release v2.0.5"), not a hand truncation.

    We run the real binary (not a mock of readglobal) because the bug is in
    Fortran sequential-read EOF semantics across two consecutive reads --
    exactly the kind of thing a unit test of parsed values, divorced from
    the actual read sequence, would not reproduce.
    """
    import shutil
    import subprocess
    import tempfile

    version = _read_version()
    exe = os.path.join(ROOT, "bin", f"run_eqdyna2d_{version}")
    fixture = os.path.join(HERE, "fixtures", "FE_Global.pre230.txt")
    if not os.path.exists(exe):
        skip("R10", "pre-2.3.0 FE_Global.txt does not crash", f"binary not built: {exe}")
        return
    if not os.path.exists(fixture):
        check("R10", "pre-2.3.0 FE_Global.txt fixture exists", False, fixture)
        return
    with tempfile.TemporaryDirectory() as d:
        shutil.copy(fixture, os.path.join(d, "FE_Global.txt"))
        r = subprocess.run([exe], cwd=d, capture_output=True, text=True, timeout=30)
    check("R10", "pre-2.3.0 FE_Global.txt: readglobal completes (exit 0)",
          r.returncode == 0,
          f"exit {r.returncode}; stdout tail: {r.stdout[-400:]}")
    check("R10", "pre-2.3.0 FE_Global.txt: no slipSense parse error",
          "could not parse slipSense" not in r.stdout,
          "the EOF-escalation regression from PR #17 is back")
    # readglobal must have returned control to the caller: the next
    # subroutine (readmodelgeometry) ran and reported its own missing file,
    # proving readglobal did not abort partway through on this fixture.
    check("R10", "pre-2.3.0 FE_Global.txt: proceeds past readglobal",
          "FE_Model_Geometry.txt is required but missing" in r.stdout,
          "did not reach readmodelgeometry; readglobal aborted or hung")


def main() -> None:
    print("EQdyna.2Dcycle convention checks (PROJECT_RULES.md)\n")
    for fn in (test_mesh_indexing,
               test_utilities_guard_index_base,
               test_nsmp_not_filtered,
               test_side_classification_uses_normal,
               test_compset_params_use_par_prefix,
               test_mesh_export_keeps_all_elements,
               test_geometry_checker,
               test_dropped_e_truncates,
               test_observed_site_constants,
               test_magnitude_constants_documented,
               test_scripts_search_araw,
               test_case_setup_run_sh,
               test_no_hardcoded_fault_counts_in_fortran,
               test_version_file_is_single_semver_line,
               test_changelog_has_body_for_version,
               test_changelog_section_extraction_returns_text,
               test_run_sh_sets_omp_threads,
               test_fig4_refuses_empty_output,
               test_stf_format_versions_match,
               test_user_guide,
               test_no_tracked_bytecode,
               test_root_whitelist,
               test_root_file_size,
               test_repo_tidy,
               test_no_unmarked_silent_fallback,
               test_status_board_exists,
               test_no_rival_status_board,
               test_feglobal_malformed_lines,
               test_fe_global_pre230_no_crash,
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
