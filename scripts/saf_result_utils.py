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

OBSERVED_EQDYNA_X_KM = np.array(
    [
        -158.576402497732,
        -69.7459994869569,
        -17.4212282973009,
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


def load_fault_blocks(case_dir: Path, fault_names: tuple[str, ...] | None = None) -> tuple[FaultBlock, ...]:
    nsmp = np.loadtxt(case_dir / "nsmp.txt", dtype=int)
    vert = np.loadtxt(case_dir / "vert.txt")
    fe_global = case_dir / "FE_Global.txt"
    c_mesh = None
    if fe_global.exists():
        with open(fe_global, "r") as f:
            first = f.readline().strip()
        if first:
            c_mesh = int(first.split()[0])
    vert_km = vert if c_mesh == 3 else vert / 1.0e3

    if nsmp.ndim != 2 or nsmp.shape[1] < 1:
        raise ValueError(f"Unexpected nsmp.txt shape: {nsmp.shape}")

    mesh_info = case_dir / "meshGeneralInfo.txt"
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
