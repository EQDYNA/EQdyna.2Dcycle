#!/usr/bin/env python3
"""Recurrence-interval and slip statistics at SSAF paleoseismic sites.

Reproduces Table 2 of Liu et al. (2022, JGR-SE, 10.1029/2021JB023420) from an
EQdyna.2Dcycle `paper.saf.A`-style run. Port of the Zenodo MATLAB script
`archive/published/zenodo_software/dunyuliu/dunyuliu-Multicycle_dynamic_SSAF_NSJF-9e645cf/
result/Figure5_Plot_Recurrene_Stats.m` (do not modify that file; it is the parity oracle).

Method (identical to the MATLAB original):
  * an event is counted at a site if fault slip at that site's node exceeds
    `--threshold` (0.5 m, per Wesnousky 2008 ~ Mw 6.0-6.25);
  * event time is the running sum of interval.txt<ic> (interseismic
    durations, years -> kyr): MATLAB ttot(1)=0, ttot(i+1) = ttot(i) +
    tinte(i+1)/1e3 -- note tinte(1) is never used, replicated verbatim;
  * recurrence intervals are successive differences of those times;
  * "mean slip" is the mean, over counted events, of the peak slip within
    the SAF node window used in the paper: slip(i, 765:1728), MATLAB
    1-based inclusive -> Python nodes [764:1728).

Site nodes are DERIVED from the mesh (vert.txt + nsmp.txt) and the paper's
`Paleo_sites_loc_inEQdyna.txt` site-coordinate file, by nearest-node search
on fault 3 (the continuous SAF/Clark trace -- every site in the MATLAB's
`plocall` table uses fault id 3). This replaces the MATLAB's hand-picked
`plocall` node indices; see LEGACY_HARDCODED below and the cross-check
printed at startup. Only valid for the paper fault-node partition
nft=(295, 178, 1769); the script refuses to guess for any other mesh.

Two known items ported verbatim without "fixing" (see MATLAB, both
irrelevant to Table 2):
  * `mom` (seismic moment) is initialized once outside the cycle loop and
    never reset per cycle, so `mag` grows monotonically and is not a
    per-event magnitude. Table 2 does not use magrec, so this bug is not
    reproduced here (no magnitude output at all).
  * the moment sum only covers fault 1 and fault 2 nodes (`jj` loops to
    nftsum(2)), never fault 3 -- again unused by Table 2.

Published Table 2 has a known error at FM (49 counted events printed vs 52
that its own script computes from its own published data). See
FM_COUNT_DISCREPANCY below; the port reports 52 and flags the delta.

Usage:
    python3 paleo_site_stats.py [case_dir] [--sites BF FM WW] [--threshold 0.5]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# Row order of archive/.../Paleo_sites_loc_inEQdyna.txt (14 rows; only the
# first 9 are named in the MATLAB's plocall / nsite comment: "nsite options
# 1,3,5,8 -> BF,FM,EL,WW" matches plocall rows 1,3,5,8 = BF,FM,EL,WW).
SITE_ORDER = ["BF", "MP", "FM", "3P", "EL", "LR", "PC", "WW", "LS"]

# plocall from Figure5_Plot_Recurrene_Stats.m: (within-fault-3 1-based index,
# fault id). Every site uses fault id 3 (the continuous SAF/Clark trace with
# 1769 nodes -- by far the largest of the 3 faults, consistent with it being
# the through-going trace all these paleo sites sit on).
LEGACY_HARDCODED = {
    "BF": 358, "MP": 33, "FM": 818, "3P": 44, "EL": 1058,
    "LR": 57, "PC": 59, "WW": 1412, "LS": 69,
}

# Fault-node partition of the paper mesh; site derivation assumes it.
PAPER_NFT = (295, 178, 1769)

# Node window over which the paper takes each event's peak slip (1-based).
MAXSLIP_WINDOW = (765, 1728)

# Default location of the paper's site-coordinate file (read-only reference
# under archive/; never modified, never copied).
_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PALEO_SITES_FILE = (
    _REPO_ROOT / "archive" / "published" / "zenodo_software" / "dunyuliu"
    / "dunyuliu-Multicycle_dynamic_SSAF_NSJF-9e645cf" / "result"
    / "Paleo_sites_loc_inEQdyna.txt"
)

# Published Table 2, Model A, for side-by-side comparison.
PAPER_TABLE2_MODEL_A = {
    #        mean(yr)  std(yr)  COV   slip(m)  slip_std(m)  counted/total
    "BF": (140.0, 86.0, 0.61, 5.7, 3.9, 107, 4000),
    "FM": (283.0, 116.0, 0.41, 8.7, 3.3, 49, 4000),   # see FM_COUNT_DISCREPANCY
    "WW": (308.0, 129.0, 0.42, 9.0, 3.2, 48, 4000),
}

# Reference MATLAB (Figure5_Plot_Recurrene_Stats.m) executed on the published
# Pangaea Model A output, full precision. This -- not the paper's rounded table
# -- is the parity oracle for this port.
MATLAB_ON_PUBLISHED_DATA = {
    #      mean(yr)     std(yr)      COV        slip(m)   slip_std(m)  n
    "BF": (140.415094, 85.965316, 0.612223, 5.731858, 3.946570, 107),
    "FM": (283.431373, 115.638273, 0.407994, 8.696537, 3.267749, 52),
    "WW": (307.553191, 128.855802, 0.418971, 8.962249, 3.222840, 48),
}

FM_COUNT_DISCREPANCY = """\
Known discrepancy in the published Table 2, at FM only.

Running the paper's own MATLAB script on the paper's own published Pangaea
output gives 52 counted events at Frazier Mountain, not the 49 printed in
Table 2. Every other FM statistic reproduces exactly (mean 283.43 -> 283,
std 115.64 -> 116, COV 0.408 -> 0.41, slip 8.697/3.268 -> 8.7/3.3), and BF
(107) and WW (48) match exactly, so this is isolated to the FM count rather
than a systematic offset.

Ruled out: node choice (swept the full fault-3 range and +/-10 nodes around
FM; no node yields 49 while keeping BF and WW exact), and threshold
sensitivity (all 52 counted events clear 0.5 m by a wide margin, min 0.74 m,
so no borderline event can be dropped by rounding).

The count therefore does not follow from the published data under the
published recipe. This port reports what the reference implementation
actually computes (52) and flags the delta, rather than tuning to 49.
"""


# --- gfortran 3-digit-exponent overflow ------------------------------------
# gfortran's E-format drops the 'E' when the exponent needs 3 digits, so
# "0.8384675E-101" is written as "0.8384675-101". np.loadtxt chokes on that.
#
# We deliberately reproduce MATLAB's reading, NOT the value gfortran intended.
# Verified against MATLAB R2024b: load() on a line "0.8384675-101 1.0 2.0"
# returns [0.8384675, 1.0, 2.0] -- its ASCII scanner takes the longest valid
# leading float and silently discards the trailing "-101". It does NOT
# reinterpret the token as scientific notation (a well-formed
# "0.5000000E-101" does read as 5e-102).
#
# So the published Table 2 / Figure 6 statistics were computed with these
# near-zero slips read as ~0.84 m. Truncating here is what makes this port
# agree with the reference implementation on the published data; interpreting
# them as E-101 would be physically right but would NOT reproduce the paper.
# Only paper-era (v2.0.2) output is affected -- current runs emit no such
# tokens -- so this rule costs nothing on new models.
_EXP_OVERFLOW = re.compile(rb"(\d)[+-]\d{3}$")


class _RepairingConverter:
    """np.loadtxt converter: truncates dropped-E tokens as MATLAB does."""

    def __init__(self) -> None:
        self.n_repaired = 0

    def __call__(self, tok: bytes) -> float:
        fixed = _EXP_OVERFLOW.sub(rb"\1", tok)
        if fixed != tok:
            self.n_repaired += 1
        return float(fixed)


def _loadtxt_repaired(path: Path, usecols) -> tuple[np.ndarray, int]:
    conv = _RepairingConverter()
    arr = np.loadtxt(path, usecols=usecols, converters={c: conv for c in usecols})
    return arr, conv.n_repaired


def _tag(path: Path) -> int:
    """Trailing icstart tag of totalop.txt<ic> / interval.txt<ic>, for ordering."""
    m = re.search(r"(\d+)$", path.name)
    return int(m.group(1)) if m else 0


def _find(case_dir: Path, pattern: str) -> list[Path]:
    """Case-dir files first, then aRawSimuData/ (run.sh moves them there).

    Segments are ordered by their numeric icstart tag, so a restarted run
    (totalop.txt1, totalop.txt1256, totalop.txt2929) stitches in time order
    rather than lexicographically.
    """
    hits = list(case_dir.glob(pattern))
    if not hits:
        hits = list((case_dir / "aRawSimuData").glob(pattern))
    return sorted(hits, key=_tag)


def _find_one(dirs: tuple[Path, ...], name: str) -> Path:
    for d in dirs:
        p = d / name
        if p.exists():
            return p
    sys.exit(f"error: {name} not found in any of {[str(d) for d in dirs]}")


def read_nft(case_dir: Path) -> tuple[int, ...]:
    # meshGeneralInfo.txt (v2.0.3+) or Mesh_general_info.txt (paper-era v2.0.2),
    # in the case dir, its aRawSimuData/, or a sibling mesh/ dir (published archive).
    names = ("meshGeneralInfo.txt", "Mesh_general_info.txt")
    dirs = (case_dir, case_dir / "aRawSimuData", case_dir.parent / "mesh")
    for d in dirs:
        for n in names:
            if (d / n).exists():
                info = d / n
                break
        else:
            continue
        break
    else:
        sys.exit(f"error: no mesh info file ({' / '.join(names)}) found near {case_dir}")
    lines = [l.split() for l in info.read_text().splitlines() if l.strip()]
    return tuple(int(v) for v in lines[1])


def derive_site_nodes(case_dir: Path, nft: tuple[int, ...],
                       sites_file: Path) -> dict[str, tuple[int, float]]:
    """Nearest-fault-3-node index (1-based within fault 3) + distance (km)
    for each site in SITE_ORDER, from vert.txt/nsmp.txt + the paper's
    site-coordinate file. Mirrors the MATLAB's x3 = vert2(nsmp2(maxftnode*2+1
    : maxftnode*2+nft(3), 1), :) extraction."""
    dirs = (case_dir, case_dir / "aRawSimuData", case_dir.parent / "mesh")
    vert = np.loadtxt(_find_one(dirs, "vert.txt")) / 1.0e3  # m -> km
    nsmp = np.loadtxt(_find_one(dirs, "nsmp.txt"), dtype=np.int64)
    if nsmp.shape[0] % len(nft):
        sys.exit(f"error: nsmp.txt has {nsmp.shape[0]} rows, not a multiple "
                  f"of {len(nft)} faults")
    maxftnode = nsmp.shape[0] // len(nft)
    fault3_ids = nsmp[maxftnode * 2: maxftnode * 2 + nft[2], 0]  # 1-based vert id
    fault3_xy = vert[fault3_ids - 1]

    if not sites_file.exists():
        sys.exit(f"error: paleo-site coordinate file not found: {sites_file}")
    pts = np.loadtxt(sites_file)
    if pts.shape[0] < len(SITE_ORDER):
        sys.exit(f"error: {sites_file} has {pts.shape[0]} rows, need at least "
                  f"{len(SITE_ORDER)} for {SITE_ORDER}")

    derived = {}
    for i, name in enumerate(SITE_ORDER):
        d = np.hypot(fault3_xy[:, 0] - pts[i, 0], fault3_xy[:, 1] - pts[i, 1])
        idx = int(np.argmin(d)) + 1  # 1-based within-fault index
        derived[name] = (idx, float(d[idx - 1]))
    return derived


def global_node(idx_1based: int, nft: tuple[int, ...]) -> int:
    """0-based global fault-node index for a fault-3 within-fault index."""
    return sum(nft[:2]) + idx_1based - 1


def load_slip(case_dir: Path, ntotnd: int) -> tuple[np.ndarray, int]:
    """(ncycles, ntotnd) slip array, cycles stacked as written by the binary."""
    files = _find(case_dir, "totalop.txt*")
    if not files:
        sys.exit(f"error: no totalop.txt* found under {case_dir}")
    blocks = []
    n_repaired_total = 0
    for f in files:
        # column 3 (1-based) is slip; read only that column to keep 800 MB files cheap
        col, n_repaired = _loadtxt_repaired(f, usecols=(2,))
        if n_repaired:
            print(f"  repaired {n_repaired} dropped-E exponent token(s) in {f.name}",
                  file=sys.stderr)
        n_repaired_total += n_repaired
        if col.size % ntotnd:
            sys.exit(f"error: {f.name} has {col.size} rows, not a multiple of ntotnd={ntotnd}")
        blocks.append(col.reshape(-1, ntotnd))
    return np.vstack(blocks), n_repaired_total


def load_times_kyr(case_dir: Path, ncycles: int) -> np.ndarray:
    """Cumulative event time in kyr; ttot[0]=0, ttot[i]=ttot[i-1]+interval[i]."""
    files = _find(case_dir, "interval.txt*")
    if not files:
        sys.exit(f"error: no interval.txt* found under {case_dir}")
    parts = []
    for f in files:
        col, n_repaired = _loadtxt_repaired(f, usecols=(0,))
        if n_repaired:
            print(f"  repaired {n_repaired} dropped-E exponent token(s) in {f.name}",
                  file=sys.stderr)
        parts.append(np.atleast_1d(col))
    iv = np.concatenate(parts)
    iv = iv[:ncycles]
    t = np.zeros(ncycles)
    t[1:] = np.cumsum(iv[1:]) / 1.0e3
    return t


def site_stats(slip: np.ndarray, t_kyr: np.ndarray, node: int,
               threshold: float, w0: int, w1: int) -> dict:
    hit = slip[:, node] > threshold
    times = t_kyr[hit]
    peak = slip[:, w0 - 1:w1].max(axis=1)[hit]
    rec = np.diff(times) * 1.0e3  # kyr -> yr
    if rec.size < 2:
        return {"n": int(hit.sum()), "mean": np.nan, "std": np.nan,
                "cov": np.nan, "slip": np.nan, "slip_std": np.nan}
    return {
        "n": int(hit.sum()),
        "mean": float(rec.mean()),
        "std": float(rec.std(ddof=1)),
        "cov": float(rec.std(ddof=1) / rec.mean()),
        "slip": float(peak.mean()),
        "slip_std": float(peak.std(ddof=1)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir", nargs="?", default=".", help="Case directory (default: cwd)")
    ap.add_argument("--sites", nargs="+", default=["BF", "FM", "WW"],
                    choices=sorted(LEGACY_HARDCODED), help="Paleoseismic sites (default: Table 2 set)")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="Minimum slip at the site to count an event, in m (default 0.5)")
    ap.add_argument("--paleo-sites-file", type=Path, default=DEFAULT_PALEO_SITES_FILE,
                    help="Site-coordinate file for node derivation (default: paper archive copy)")
    ap.add_argument("--csv", default=None, help="Also write the table to this CSV path")
    args = ap.parse_args()

    case_dir = Path(args.case_dir).resolve()
    nft = read_nft(case_dir)
    if tuple(nft) != PAPER_NFT:
        sys.exit(f"error: site derivation assumes the paper mesh nft={PAPER_NFT}, "
                 f"but this case has nft={nft}. Re-derive site nodes before using them.")
    ntotnd = sum(nft)

    derived = derive_site_nodes(case_dir, nft, args.paleo_sites_file)
    print("site node derivation (fault 3, 1-based) -- derived vs legacy MATLAB plocall:")
    print(f"  {'site':<5} {'derived':>8} {'dist(km)':>9} {'legacy':>7} {'match':>6}")
    for s in SITE_ORDER:
        idx, dist = derived[s]
        legacy = LEGACY_HARDCODED[s]
        match = "yes" if idx == legacy else f"no(+{idx - legacy})"
        print(f"  {s:<5} {idx:>8d} {dist:>9.3f} {legacy:>7d} {match:>6}")
    print()

    slip, n_repaired = load_slip(case_dir, ntotnd)
    ncycles = slip.shape[0]
    t_kyr = load_times_kyr(case_dir, ncycles)
    print(f"case      : {case_dir}")
    print(f"cycles    : {ncycles}   simulated: {t_kyr[-1]:.2f} kyr   "
          f"(paper Model A: 4000 cycles / 15 kyr)")
    print(f"threshold : {args.threshold} m slip at site")
    print(f"dropped-E tokens truncated (MATLAB rule): {n_repaired}\n")

    hdr = (f"{'site':<5} {'mean rec (yr)':>14} {'COV':>6} {'slip (m)':>9} "
           f"{'sd slip (m)':>12} {'counted/total':>14}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for s in args.sites:
        idx, _ = derived[s]
        node = global_node(idx, nft)
        st = site_stats(slip, t_kyr, node, args.threshold, *MAXSLIP_WINDOW)
        print(f"{s:<5} {st['mean']:>8.0f}(+/-{st['std']:.0f}) {st['cov']:>6.2f} "
              f"{st['slip']:>9.1f} {st['slip_std']:>12.1f} {st['n']:>9d}/{ncycles}")
        rows.append((s, st))
        if s in PAPER_TABLE2_MODEL_A:
            m, sd, cov, sl, slsd, n, tot = PAPER_TABLE2_MODEL_A[s]
            print(f"{'  paper':<5} {m:>8.0f}(+/-{sd:.0f}) {cov:>6.2f} "
                  f"{sl:>9.1f} {slsd:>12.1f} {n:>9d}/{tot}")

    # Flag the known Table 2 error rather than letting it read as a port bug.
    if "FM" in args.sites:
        print("\n" + "-" * len(hdr))
        print(FM_COUNT_DISCREPANCY.rstrip())

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("site,mean_recurrence_yr,std_recurrence_yr,cov,mean_slip_m,std_slip_m,n_events,n_cycles\n")
            for s, st in rows:
                f.write(f"{s},{st['mean']:.2f},{st['std']:.2f},{st['cov']:.4f},"
                        f"{st['slip']:.3f},{st['slip_std']:.3f},{st['n']},{ncycles}\n")
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
