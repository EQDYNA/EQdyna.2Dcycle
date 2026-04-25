#! /usr/bin/env python3
"""Quick smoke test: build + create case + run 1 cycle + check totalop.txt1 exists.

Use as a fast sanity check before running the full test suite.
"""
import glob, os, sys, subprocess, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
os.chdir(ROOT)

print("Building...")
if subprocess.call(['./install.sh', '-m', 'ubuntu']) != 0:
    print("build failed"); sys.exit(1)

os.environ['EQDYNA2DCYCLEROOT'] = ROOT
os.environ['PATH'] = f"{ROOT}/bin:{ROOT}/scripts:{os.environ['PATH']}"

case = 'work/smoke_case'
subprocess.call(['rm', '-rf', case])
print(f"🏗️  Creating case at {case}...")
if subprocess.call(['python3', 'scripts/create.newcase',
                    '--work_dir', case, '--compset', 'subei.gmsh.lite', '--force']) != 0:
    print("❌ create.newcase failed"); sys.exit(1)

os.chdir(case)

# Force a single cycle: smoke just verifies the binary runs end-to-end.
import re
with open('user_defined_params.py') as f: src = f.read()
src = re.sub(r'par\.icstart\s*,\s*par\.icend\s*=.*',
             'par.icstart, par.icend = 1, 1', src)
with open('user_defined_params.py', 'w') as f: f.write(src)

print("🌐 meshgen + case.setup...")
if subprocess.call(['python3', 'meshgen.py']) != 0: sys.exit(1)
if subprocess.call(['python3', 'case.setup']) != 0: sys.exit(1)

print("run.sh (1 cycle)...")
if subprocess.call(['bash', 'run.sh']) != 0:
    print("run.sh failed to launch"); sys.exit(1)

# run.sh backgrounds the simulation via nohup; poll the log until it prints
# "Job finished" or we time out.
deadline = time.time() + 600
while time.time() < deadline:
    logs = sorted(glob.glob('run_*.log'))
    if logs:
        with open(logs[-1]) as f:
            if 'Job finished' in f.read():
                break
    time.sleep(5)
else:
    print("timed out waiting for run.sh to finish (>600s)"); sys.exit(1)

out = 'aRawSimuData/totalop.txt1'
if not os.path.exists(out) or os.path.getsize(out) == 0:
    print(f"{out} missing or empty"); sys.exit(1)

print(f"smoke test passed ({out} = {os.path.getsize(out)} bytes)")
