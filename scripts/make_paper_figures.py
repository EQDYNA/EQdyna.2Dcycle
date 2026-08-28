#!/usr/bin/env python3
"""Generate the full Liu et al. (2022) figure set for any EQdyna.2Dcycle case.

One command that produces, for a new model, every figure and table the paper
presents -- so a new run can be read against the published results without
remembering which script does what.

Stages (each maps to a published figure or table):

    catalog     CATALOG=1 plotRuptureDynamics  -> aPlots/catalog.csv
                per-event magnitude, moment, nucleation, duration, peak slip
    rupture     plotRuptureDynamics             -> aPlots/cRuptureDynamics*.png
                per-cycle 4-panel shear/normal/slip/rupture-time
    figure3     plot_saf_figure3.py             -> long-term slip rates vs UCERF3
    figure4     plot_event_slips_overtime_fig4  -> Fig 4/7/8 slip-distribution stacks
    figure5     paleo_site_stats.py             -> Fig 5 / Table 2 recurrence stats
    figure6     plot_saf_figure6.py             -> moment release + magnitude-frequency
    figure9     plot_saf_figure9.py             -> characteristic-event slip distributions
    analysis    analyze_catalog.py              -> b-value, MFD, magnitude-vs-cycle

Stages whose script is not present yet are reported as SKIP, not as failure,
so this driver stays usable while the remaining MATLAB ports land.

Usage:
    python3 make_paper_figures.py [case_dir] [--only figure4 figure5] [--skip rupture]
                                  [--min-magnitude 6.5] [--force]

Notes:
  * `rupture` is by far the most expensive stage (one figure per event above
    --min-magnitude). --skip rupture when iterating on the summary figures.
  * Existing figures are left alone unless --force is given.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Stage name -> (script filename, extra argv, passes case_dir as argv, env overrides)
# Scripts are looked up next to this file first, then in the case dir (cases get
# their own copies via create.newcase), then on PATH.
STAGES = [
    ("catalog",  "plotRuptureDynamics",            [], False, {"CATALOG": "1"}),
    ("figure3",  "plot_saf_figure3.py",            [], True,  {}),
    ("figure4",  "plot_event_slips_overtime_fig4.py", [], True, {}),
    ("figure5",  "paleo_site_stats.py",            [], True,  {}),
    ("figure6",  "plot_saf_figure6.py",            [], True,  {}),
    ("figure9",  "plot_saf_figure9.py",            [], True,  {}),
    ("analysis", "analyze_catalog.py",             [], True,  {}),
    # rupture last: it is the expensive one, so the cheap summary figures
    # are already on disk if the user interrupts.
    ("rupture",  "plotRuptureDynamics",            [], False, {}),
]

SCRIPT_DIR = Path(__file__).resolve().parent


def find_script(name: str, case_dir: Path) -> Path | None:
    for cand in (SCRIPT_DIR / name, case_dir / name):
        if cand.exists():
            return cand
    found = shutil.which(name)
    return Path(found) if found else None


def run_stage(stage: str, script: Path, extra: list[str], pass_case: bool,
              env_over: dict, case_dir: Path, log_dir: Path) -> tuple[str, float, str]:
    argv = [sys.executable, str(script)]
    if pass_case:
        argv.append(str(case_dir))
    argv += extra

    env = dict(os.environ)
    env.update(env_over)

    log_path = log_dir / f"{stage}.log"
    t0 = time.time()
    with open(log_path, "w") as log:
        proc = subprocess.run(argv, cwd=case_dir, env=env, stdout=log,
                              stderr=subprocess.STDOUT)
    dt = time.time() - t0
    status = "OK" if proc.returncode == 0 else f"FAIL(rc={proc.returncode})"
    return status, dt, str(log_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir", nargs="?", default=".", help="Case directory (default: cwd)")
    ap.add_argument("--only", nargs="+", metavar="STAGE",
                    help="Run only these stages")
    ap.add_argument("--skip", nargs="+", metavar="STAGE", default=[],
                    help="Skip these stages (e.g. --skip rupture)")
    ap.add_argument("--min-magnitude", type=float, default=None,
                    help="Magnitude gate for per-cycle rupture figures "
                         "(default: plotRuptureDynamics' own default of 6.5)")
    ap.add_argument("--force", action="store_true",
                    help="Re-render figures that already exist")
    ap.add_argument("--list", action="store_true", help="List stages and exit")
    args = ap.parse_args()

    if args.list:
        for name, script, _, _, _ in STAGES:
            print(f"  {name:<9} {script}")
        return

    case_dir = Path(args.case_dir).resolve()
    if not case_dir.is_dir():
        sys.exit(f"error: {case_dir} is not a directory")

    plots_dir = case_dir / "aPlots"
    plots_dir.mkdir(exist_ok=True)
    log_dir = plots_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    selected = [s for s in STAGES
                if (args.only is None or s[0] in args.only) and s[0] not in args.skip]
    if args.only:
        unknown = set(args.only) - {s[0] for s in STAGES}
        if unknown:
            sys.exit(f"error: unknown stage(s): {', '.join(sorted(unknown))}")

    print(f"case   : {case_dir}")
    print(f"output : {plots_dir}")
    print(f"stages : {', '.join(s[0] for s in selected)}\n")

    results = []
    for stage, script_name, extra, pass_case, env_over in selected:
        script = find_script(script_name, case_dir)
        if script is None:
            print(f"  {stage:<9} SKIP  ({script_name} not found -- port not landed yet?)")
            results.append((stage, "SKIP", 0.0))
            continue

        env_over = dict(env_over)
        if stage == "rupture":
            if args.min_magnitude is not None:
                env_over["MIN_PLOT_MAGNITUDE"] = str(args.min_magnitude)
            if args.force:
                env_over["FORCE_REPLOT"] = "1"

        print(f"  {stage:<9} running {script.name} ...", end="", flush=True)
        status, dt, log_path = run_stage(stage, script, extra, pass_case,
                                         env_over, case_dir, log_dir)
        print(f"\r  {stage:<9} {status:<12} {dt:7.1f}s   log: {log_path}")
        results.append((stage, status, dt))

    print()
    ok = sum(1 for _, s, _ in results if s == "OK")
    skipped = sum(1 for _, s, _ in results if s == "SKIP")
    failed = [n for n, s, _ in results if s.startswith("FAIL")]
    print(f"{ok} ok, {skipped} skipped, {len(failed)} failed"
          + (f" ({', '.join(failed)})" if failed else ""))
    if failed:
        print("Failed stages left their output in aPlots/logs/<stage>.log")
        sys.exit(1)


if __name__ == "__main__":
    main()
