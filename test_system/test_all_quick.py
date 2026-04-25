#! /usr/bin/env python3
"""
EQdyna.2Dcycle quick smoke test.

Runs every compset (excluding paper.saf.A) for icend=10 and post-processes
the result with catalog + b-value analysis + rupture-dynamics + Figure 4.

Use this for fast end-to-end verification (build -> mesh -> simulate -> plot).
For full convergence runs use test_all.py.
"""
import glob
import os
import re
import sys
import time

QUICK_ICEND = 10
WAIT_TIMEOUT_SEC = 1800
POLL_SEC = 5

COMPSETS = [
    "saf.gmsh.lite",
    "subei.gmsh.lite",
    "gulang.gmsh.lite",
    "test.saf",
    "test.subei",
    "test.gulang",
]

print("EQdyna.2Dcycle quick smoke test")
print("=" * 50)

if os.path.basename(os.getcwd()) in ("test", "test_system"):
    os.chdir("..")
    print(f"Changed to project root: {os.getcwd()}")

print("Cleaning previous quick-test results...")
os.system("rm -rf work/test_quick_results")

print("Building EQdyna.2Dcycle...")
if os.system("./install.sh -m ubuntu") != 0:
    sys.exit("Build failed.")

os.environ["EQDYNA2DCYCLEROOT"] = os.getcwd()
os.environ["PATH"] = f"{os.getcwd()}/bin:{os.getcwd()}/scripts:{os.environ['PATH']}"

os.makedirs("work/test_quick_results", exist_ok=True)
os.chdir("work/test_quick_results")
results_root = os.getcwd()

start_time = time.time()


def patch_icend(params_path: str, icend: int) -> None:
    with open(params_path) as f:
        src = f.read()
    new = re.sub(
        r"par\.icstart\s*,\s*par\.icend\s*=.*",
        f"par.icstart, par.icend = 1, {icend}",
        src,
    )
    if new == src:
        new = src + f"\npar.icstart, par.icend = 1, {icend}\n"
    with open(params_path, "w") as f:
        f.write(new)


def wait_for_run(case_dir: str) -> bool:
    deadline = time.time() + WAIT_TIMEOUT_SEC
    while time.time() < deadline:
        logs = sorted(glob.glob(os.path.join(case_dir, "run_*.log")))
        if logs:
            with open(logs[-1]) as f:
                if "Job finished" in f.read():
                    return True
        time.sleep(POLL_SEC)
    return False


def post_process(case_dir: str) -> None:
    cwd = os.getcwd()
    os.chdir(case_dir)
    print("  catalog...")
    os.system("CATALOG=1 python3 ../../../scripts/plotRuptureDynamics > /dev/null 2>&1")
    print("  b-value analysis...")
    os.system("python3 ../../../scripts/analyze_catalog.py . --mmax 7.0")
    print("  rupture-dynamics plots...")
    os.system("MIN_PLOT_MAGNITUDE=0 python3 ../../../scripts/plotRuptureDynamics > /dev/null 2>&1")
    print("  Figure 4 slip-over-time...")
    os.system("python3 ../../../scripts/plot_event_slips_overtime_fig4.py . --threshold 0 --duration 1.0")
    os.chdir(cwd)


def run_case(compset: str) -> bool:
    print(f"\n--- {compset} ---")
    if os.system(f"python3 ../../scripts/create.newcase --work_dir {compset} --compset {compset} --force") != 0:
        print(f"  create.newcase failed for {compset}")
        return False
    patch_icend(os.path.join(compset, "user_defined_params.py"), QUICK_ICEND)
    os.chdir(compset)
    case_dir = os.getcwd()
    if os.system("python3 case.setup") != 0:
        print(f"  case.setup failed for {compset}")
        os.chdir(results_root)
        return False
    if os.system("bash run.sh") != 0:
        print(f"  run.sh failed to launch for {compset}")
        os.chdir(results_root)
        return False
    os.chdir(results_root)
    if not wait_for_run(case_dir):
        print(f"  TIMEOUT waiting for {compset} (>{WAIT_TIMEOUT_SEC}s)")
        return False
    post_process(case_dir)
    print(f"  {compset} OK")
    return True


print(f"\nRunning {len(COMPSETS)} compset(s) at icend={QUICK_ICEND}:")
ok = sum(1 for c in COMPSETS if run_case(c))

elapsed = time.time() - start_time
print(f"\nSummary: {ok}/{len(COMPSETS)} passed in {elapsed:.1f}s")
sys.exit(0 if ok == len(COMPSETS) else 1)
