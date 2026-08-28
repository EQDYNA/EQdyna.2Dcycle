#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


FAULT_NAMES = ("sjfn", "sjfs", "ssaf")

OBSERVED_RATE_ROWS = np.array(
    [
        [31.0, 37.0, 34.0],  # Carrizo
        [20.0, 34.0, 27.0],  # Mojave S
        [5.0, 20.0, 13.0],   # San Bernardino S
        [2.0, 10.0, 6.0],    # NSJF branch 1
        [6.0, 11.0, 8.0],    # NSJF branch 2
    ]
)

# Along-strike EQdyna x positions (km) of the three SSAF observed-rate sites
# (Carrizo, Mojave S, San Bernardino S) from Dawson & Weldon (2003).
#
# Derived by running the published MATLAB chain, not by hand: the lon/lat in
# observed_sliprates.m through deg2utm -> convert(x0=3.7e5, y0=3.8e6) ->
# rotate(theta=40). Verified against MATLAB R2024b, which returns
#   Carrizo          -164.315662175388
#   Mojave S           21.621033567916
#   San Bernardino S  129.481611408865
#
# A previous unverified set of values here (-158.576, -69.746, -17.421) put the
# latter two roughly 91 km and 147 km too far NW, which misassigned about a
# third of the SSAF nodes to the wrong observed-rate band. Do not change these
# without re-running the MATLAB chain.
OBSERVED_EQDYNA_X_KM = np.array(
    [
        -164.315662175388,
        21.621033567916,
        129.481611408865,
    ]
)


@dataclass(frozen=True)
class FaultBlock:
    name: str
    x_km: np.ndarray
    y_km: np.ndarray
    index_start: int
    index_stop: int

    @property
    def count(self) -> int:
        return self.index_stop - self.index_start


@dataclass(frozen=True)
class SafCaseData:
    case_dir: Path
    source_dir: Path
    cycle_tags: tuple[int, ...]
    cycle_start: int
    cycle_end: int
    intervals_kyr: np.ndarray
    event_times_kyr: np.ndarray
    slips_m: np.ndarray
    fault_blocks: tuple[FaultBlock, ...]

    @property
    def total_fault_nodes(self) -> int:
        return self.slips_m.shape[1]

    @property
    def total_events(self) -> int:
        return self.slips_m.shape[0]


def _numeric_suffix(path: Path, stem: str) -> int | None:
    match = re.fullmatch(rf"{re.escape(stem)}(\d+)", path.name)
    if not match:
        return None
    return int(match.group(1))


def discover_cycle_tags(case_dir: Path) -> tuple[Path, tuple[int, ...]]:
    search_dirs = [case_dir, case_dir / "aRawSimuData"]
    best_dir = None
    best_tags: tuple[int, ...] = ()

    for root in search_dirs:
        if not root.exists():
            continue
        tags = []
        for path in root.glob("totalop.txt*"):
            tag = _numeric_suffix(path, "totalop.txt")
            if tag is None:
                continue
            if (root / f"cyclelog.txt{tag}").exists() and (root / f"interval.txt{tag}").exists():
                tags.append(tag)
        tags_tuple = tuple(sorted(tags))
        if len(tags_tuple) > len(best_tags):
            best_dir = root
            best_tags = tags_tuple

    if best_dir is None or not best_tags:
        raise FileNotFoundError(f"No cyclelog/interval/totalop triplets found under {case_dir}")

    return best_dir, best_tags


def load_fault_blocks(
    case_dir: Path,
    fault_names: tuple[str, ...] | None = None,
    mesh_dir: Path | None = None,
) -> tuple[FaultBlock, ...]:
    md = mesh_dir if mesh_dir is not None else case_dir
    nsmp = np.loadtxt(md / "nsmp.txt", dtype=int)
    vert = np.loadtxt(md / "vert.txt")
    fe_global = md / "FE_Global.txt"
    c_mesh = None
    if fe_global.exists():
        with open(fe_global, "r") as f:
            first = f.readline().strip()
        if first:
            c_mesh = int(first.split()[0])
    vert_km = vert if c_mesh == 3 else vert / 1.0e3

    if nsmp.ndim != 2 or nsmp.shape[1] < 1:
        raise ValueError(f"Unexpected nsmp.txt shape: {nsmp.shape}")

    mesh_info = md / "meshGeneralInfo.txt"
    if not mesh_info.exists():
        mesh_info = md / "Mesh_general_info.txt"
    declared_counts: list[int] | None = None
    if mesh_info.exists():
        with open(mesh_info, "r") as f:
            next(f)
            declared_counts = list(map(int, f.readline().split()))

    if fault_names is None:
        ntotft = len(declared_counts) if declared_counts else len(FAULT_NAMES)
        fault_names = FAULT_NAMES if ntotft == len(FAULT_NAMES) else tuple(f"ft{i+1}" for i in range(ntotft))

    if declared_counts is not None and len(declared_counts) != len(fault_names):
        declared_counts = None
    max_fault_nodes = nsmp.shape[0] // len(fault_names)

    blocks = []
    running_start = 0

    for block_id, fault_name in enumerate(fault_names):
        block_rows = nsmp[block_id * max_fault_nodes:(block_id + 1) * max_fault_nodes, 0]
        if declared_counts is not None:
            node_ids = block_rows[:declared_counts[block_id]]
        else:
            # C_mesh=2 uses 1-based ids with zero padding; C_mesh=3 uses 0-based
            # ids with zero padding. If any negative ids exist, something is wrong.
            if np.any(block_rows < 0):
                raise ValueError("Unexpected negative node ids in nsmp.txt")
            # Fallback assumes the active rows are contiguous from the top.
            last_nonzero = np.max(np.flatnonzero((block_rows != 0) | (np.arange(block_rows.size) == 0))) + 1
            node_ids = block_rows[:last_nonzero]
            if np.min(node_ids) >= 1:
                node_ids = node_ids - 1
        xy = vert_km[node_ids]
        index_stop = running_start + xy.shape[0]
        blocks.append(
            FaultBlock(
                name=fault_name,
                x_km=xy[:, 0],
                y_km=xy[:, 1],
                index_start=running_start,
                index_stop=index_stop,
            )
        )
        running_start = index_stop

    return tuple(blocks)


def load_saf_case(case_dir: str | Path, cycle_tags: list[int] | tuple[int, ...] | None = None) -> SafCaseData:
    case_dir = Path(case_dir)
    fault_blocks = load_fault_blocks(case_dir)
    total_fault_nodes = sum(block.count for block in fault_blocks)

    if cycle_tags:
        source_dir = case_dir if (case_dir / f"totalop.txt{cycle_tags[0]}").exists() else case_dir / "aRawSimuData"
        tags = tuple(int(tag) for tag in cycle_tags)
    else:
        source_dir, tags = discover_cycle_tags(case_dir)

    cycle_start = None
    cycle_end = None
    intervals_parts = []
    slip_parts = []

    for tag in tags:
        cyclelog = np.atleast_1d(np.loadtxt(source_dir / f"cyclelog.txt{tag}")).astype(int)
        if cyclelog.size < 2:
            raise ValueError(f"cyclelog.txt{tag} should contain at least two integers")

        if cycle_start is None:
            cycle_start = int(cyclelog[0])
        cycle_end = int(cyclelog[1])

        intervals = np.atleast_1d(np.loadtxt(source_dir / f"interval.txt{tag}")).astype(float)
        slips = np.atleast_1d(
            np.loadtxt(source_dir / f"totalop.txt{tag}", usecols=(2,), dtype=float)
        )
        usable_rows = (slips.size // total_fault_nodes) * total_fault_nodes
        if usable_rows == 0:
            raise ValueError(f"totalop.txt{tag} does not contain a complete event block")
        slips = slips[:usable_rows]
        slip_parts.append(slips.reshape((-1, total_fault_nodes)))
        intervals_parts.append(intervals.reshape(-1))

    assert cycle_start is not None and cycle_end is not None

    slips_m = np.vstack(slip_parts)
    intervals_kyr = np.concatenate(intervals_parts) / 1.0e3
    event_times_kyr = np.zeros(intervals_kyr.shape[0], dtype=float)
    if intervals_kyr.shape[0] > 1:
        event_times_kyr[1:] = np.cumsum(intervals_kyr[1:])

    expected_events = cycle_end - cycle_start + 1
    available_events = min(expected_events, slips_m.shape[0], event_times_kyr.shape[0])
    slips_m = slips_m[:available_events]
    event_times_kyr = event_times_kyr[:available_events]
    cycle_end = cycle_start + available_events - 1

    return SafCaseData(
        case_dir=case_dir,
        source_dir=source_dir,
        cycle_tags=tags,
        cycle_start=cycle_start,
        cycle_end=cycle_end,
        intervals_kyr=intervals_kyr,
        event_times_kyr=event_times_kyr,
        slips_m=slips_m,
        fault_blocks=fault_blocks,
    )


def geologic_rate_profile(block: FaultBlock) -> np.ndarray:
    rates = np.zeros((block.count, 3), dtype=float)
    if block.name == "ssaf":
        midpoint_1 = (OBSERVED_EQDYNA_X_KM[1] - OBSERVED_EQDYNA_X_KM[0]) / 2.0 + OBSERVED_EQDYNA_X_KM[0]
        midpoint_2 = (OBSERVED_EQDYNA_X_KM[2] - OBSERVED_EQDYNA_X_KM[1]) / 2.0 + OBSERVED_EQDYNA_X_KM[1]
        left = block.x_km < midpoint_1
        middle = (block.x_km >= midpoint_1) & (block.x_km < midpoint_2)
        right = block.x_km >= midpoint_2
        rates[left] = OBSERVED_RATE_ROWS[0]
        rates[middle] = OBSERVED_RATE_ROWS[1]
        rates[right] = OBSERVED_RATE_ROWS[2]
    elif block.name == "sjfn":
        rates[:] = OBSERVED_RATE_ROWS[3]
    elif block.name == "sjfs":
        rates[:] = OBSERVED_RATE_ROWS[4]
    else:
        raise ValueError(f"Unknown fault block name: {block.name}")
    return rates


def modeled_slip_rate_mm_per_yr(case_data: SafCaseData) -> np.ndarray:
    if case_data.total_events < 2:
        raise ValueError("Need at least two events to compute long-term slip rate")
    tspan_kyr = case_data.event_times_kyr[-1] - case_data.event_times_kyr[0]
    if tspan_kyr <= 0:
        raise ValueError(f"Non-positive time span: {tspan_kyr}")
    cumulative_slip_m = np.sum(case_data.slips_m, axis=0)
    return cumulative_slip_m / tspan_kyr


def select_time_window(case_data: SafCaseData, tstart_kyr: float, tend_kyr: float) -> tuple[np.ndarray, np.ndarray]:
    mask = (case_data.event_times_kyr >= tstart_kyr) & (case_data.event_times_kyr <= tend_kyr)
    indices = np.flatnonzero(mask)
    cycle_ids = case_data.cycle_start + indices
    return indices, cycle_ids


# --- mesh-dir resolution for archive-layout cases --------------------------
# Published Pangaea work_vis*/ case dirs keep totalop.txt* etc alongside the
# cycle logs but the mesh (nsmp.txt, vert.txt, meshGeneralInfo.txt) lives in
# a SIBLING mesh/ directory, not inside case_dir. Current-format cases keep
# the mesh files in case_dir (or its aRawSimuData/). Try both.
def resolve_mesh_dir(case_dir: Path) -> Path:
    for d in (case_dir, case_dir / "aRawSimuData", case_dir.parent / "mesh"):
        if (d / "nsmp.txt").exists() and (d / "vert.txt").exists():
            return d
    raise FileNotFoundError(f"No nsmp.txt/vert.txt found for case near {case_dir}")


# --- gfortran 3-digit-exponent overflow repair ------------------------------
# gfortran's E-format drops the 'E' when the exponent needs 3 digits, e.g.
# "0.8384675E-101" is written as "0.8384675-101". np.loadtxt chokes on that.
#
# IMPORTANT, and a DIVERGENCE from scripts/paleo_site_stats.py's repair: this
# is NOT reinterpreted as scientific notation ("...E-101"). Verified against
# MATLAB's own `load()` on a synthetic single-token file
# ("0.8384675-101  1.0  2.0" -> x = [0.8384675 1.0 2.0], NOT
# [0.8384675e-101 1.0 2.0]): MATLAB's ASCII-matrix scanner reads the longest
# valid leading float from each whitespace token and silently discards the
# unparsed trailing "-101" (no extra column is created, no error is raised).
# The parity oracle for Figure6/Figure9 is that behavior -- MATLAB's `load`
# on the real archive data confirmed the truncated value (e.g. slip(1,:) sum
# for work_vis7_fs0.5/totalop.txt1 event 1 comes out to 20.4948 under
# truncation vs 3.3657 under exponent-reinsertion; MATLAB's own `load`
# reports 20.4948). paleo_site_stats.py's exponent-reinsertion repair was
# never exercised on a divergence-sensitive column/row (its analysed nodes
# happened not to hit an affected token), so its Table 2 numbers still came
# out matching MATLAB's -- but that repair does NOT reproduce MATLAB's
# actual parsing rule and would silently diverge on data that does hit one.
# Not fixed here (out of scope / different file ownership); see the port
# report.
_EXP_OVERFLOW = re.compile(rb"(\d)([+-]\d{3})$")


class _RepairingConverter:
    """np.loadtxt converter: truncates dropped-E tokens to their leading
    valid mantissa (matching MATLAB `load()`'s scan-and-stop behavior),
    counts how many tokens were affected."""

    def __init__(self) -> None:
        self.n_repaired = 0

    def __call__(self, tok: bytes) -> float:
        fixed = _EXP_OVERFLOW.sub(rb"\1", tok)
        if fixed != tok:
            self.n_repaired += 1
        return float(fixed)


def loadtxt_repaired(path: Path, usecols) -> tuple[np.ndarray, int]:
    conv = _RepairingConverter()
    arr = np.loadtxt(path, usecols=usecols, converters={c: conv for c in usecols})
    return arr, conv.n_repaired


# --- shared moment / magnitude kernel (Figure6 & Figure9) ------------------
# Ports the mom/len accumulation loops that are byte-identical across
# Figure6_Plot_Magnitude_Frequency.m and Figure9_Special_Events.m. Per fault
# block (block boundaries are NOT bridged -- each block's end nodes only get
# a one-sided half-segment weight), every node gets a trapezoidal
# along-strike weight: half the length of each adjacent segment. mom sums
# slip*weight*rig (weight still in km at this point); len sums weight only
# over nodes with positive slip; the final *1e3 converts the accumulated
# km-weights to metres in one shot (equivalent to converting every km
# distance to metres inside the sum, since the sum is linear).
def fault_node_weights_km(fault_blocks: tuple[FaultBlock, ...]) -> np.ndarray:
    total = sum(b.count for b in fault_blocks)
    w = np.zeros(total, dtype=float)
    for b in fault_blocks:
        n = b.count
        if n == 0:
            continue
        if n == 1:
            continue  # w[block] stays 0, matching the MATLAB loop's implicit
                      # behaviour for a degenerate single-node block (no jj-1
                      # or jj+1 branch ever fires for the sole node).
        seglen = np.hypot(np.diff(b.x_km), np.diff(b.y_km))
        bw = np.empty(n, dtype=float)
        bw[0] = 0.5 * seglen[0]
        bw[-1] = 0.5 * seglen[-1]
        bw[1:-1] = 0.5 * (seglen[:-1] + seglen[1:])
        w[b.index_start:b.index_stop] = bw
    return w


def moment_and_magnitude(
    slips_m: np.ndarray,
    weights_km: np.ndarray,
    rig: float,
    lockdepth_m: float,
    alpha: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized port of the mom/mag accumulation in Figure6/Figure9.

    slips_m: (nevents, ntotnd). Returns (mom [N*m], mag [Mw]) each (nevents,).
    """
    mom_raw = rig * (slips_m @ weights_km)
    len_km = (slips_m > 0.0) @ weights_km
    mom = alpha * mom_raw * np.minimum(len_km * 1.0e3, lockdepth_m) * 1.0e3
    mag = (2.0 / 3.0) * np.log10(mom) - 6.07
    return mom, mag


# --- full-column cycle loader (ss, ns, slip, rupt) for Figure6 / Figure9 ---
@dataclass(frozen=True)
class SafFullCaseData:
    case_dir: Path
    source_dir: Path
    mesh_dir: Path
    cycle_tags: tuple[int, ...]
    cycle_start: int
    cycle_end: int
    event_times_kyr: np.ndarray
    shear_stress: np.ndarray   # (nevents, ntotnd), totalop col 1 (1-based)
    normal_stress: np.ndarray  # totalop col 2
    slip_m: np.ndarray         # totalop col 3
    rupture_time_s: np.ndarray  # totalop col 5
    fault_blocks: tuple[FaultBlock, ...]

    @property
    def total_fault_nodes(self) -> int:
        return self.slip_m.shape[1]

    @property
    def total_events(self) -> int:
        return self.slip_m.shape[0]


def load_saf_full_case(
    case_dir: str | Path,
    cycle_tags: list[int] | tuple[int, ...] | None = None,
) -> SafFullCaseData:
    """Loads ss/ns/slip/rupt for every cycle, mirroring the load loop shared
    by Figure6_Plot_Magnitude_Frequency.m and Figure9_Special_Events.m:
    concatenate totalop.txt<tag>/interval.txt<tag> segments in tag order,
    cycle_start from the first segment's cyclelog, cycle_end from the last."""
    case_dir = Path(case_dir)
    mesh_dir = resolve_mesh_dir(case_dir)
    fault_blocks = load_fault_blocks(case_dir, mesh_dir=mesh_dir)
    total_fault_nodes = sum(block.count for block in fault_blocks)

    if cycle_tags:
        source_dir = case_dir if (case_dir / f"totalop.txt{cycle_tags[0]}").exists() else case_dir / "aRawSimuData"
        tags = tuple(int(t) for t in cycle_tags)
    else:
        source_dir, tags = discover_cycle_tags(case_dir)

    cycle_start = None
    cycle_end = None
    intervals_parts = []
    ss_parts, ns_parts, slip_parts, rupt_parts = [], [], [], []
    n_repaired_total = 0

    for tag in tags:
        cyclelog = np.atleast_1d(np.loadtxt(source_dir / f"cyclelog.txt{tag}")).astype(int)
        if cyclelog.size < 2:
            raise ValueError(f"cyclelog.txt{tag} should contain at least two integers")
        if cycle_start is None:
            cycle_start = int(cyclelog[0])
        cycle_end = int(cyclelog[1])

        intervals, n_rep = loadtxt_repaired(source_dir / f"interval.txt{tag}", usecols=(0,))
        n_repaired_total += n_rep
        intervals_parts.append(np.atleast_1d(intervals))

        cols, n_rep = loadtxt_repaired(source_dir / f"totalop.txt{tag}", usecols=(0, 1, 2, 4))
        n_repaired_total += n_rep
        usable_rows = (cols.shape[0] // total_fault_nodes) * total_fault_nodes
        if usable_rows == 0:
            raise ValueError(f"totalop.txt{tag} does not contain a complete event block")
        cols = cols[:usable_rows].reshape(-1, total_fault_nodes, 4)
        ss_parts.append(cols[:, :, 0])
        ns_parts.append(cols[:, :, 1])
        slip_parts.append(cols[:, :, 2])
        rupt_parts.append(cols[:, :, 3])

    assert cycle_start is not None and cycle_end is not None

    ss = np.vstack(ss_parts)
    ns = np.vstack(ns_parts)
    slip_m = np.vstack(slip_parts)
    rupt = np.vstack(rupt_parts)
    intervals_kyr = np.concatenate(intervals_parts) / 1.0e3
    event_times_kyr = np.zeros(intervals_kyr.shape[0], dtype=float)
    if intervals_kyr.shape[0] > 1:
        event_times_kyr[1:] = np.cumsum(intervals_kyr[1:])

    expected_events = cycle_end - cycle_start + 1
    available_events = min(expected_events, slip_m.shape[0], event_times_kyr.shape[0])
    ss = ss[:available_events]
    ns = ns[:available_events]
    slip_m = slip_m[:available_events]
    rupt = rupt[:available_events]
    event_times_kyr = event_times_kyr[:available_events]
    cycle_end = cycle_start + available_events - 1

    if n_repaired_total:
        import sys
        print(f"  repaired {n_repaired_total} dropped-E exponent token(s) under {source_dir}",
              file=sys.stderr)

    return SafFullCaseData(
        case_dir=case_dir,
        source_dir=source_dir,
        mesh_dir=mesh_dir,
        cycle_tags=tags,
        cycle_start=cycle_start,
        cycle_end=cycle_end,
        event_times_kyr=event_times_kyr,
        shear_stress=ss,
        normal_stress=ns,
        slip_m=slip_m,
        rupture_time_s=rupt,
        fault_blocks=fault_blocks,
    )
