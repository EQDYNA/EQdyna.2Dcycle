#! /usr/bin/env python3
"""Quick smoke test: build + create case + run 1 cycle + check totalop.txt1 exists.

Use as a fast sanity check before running the full test suite.
"""
import os, sys, subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
os.chdir(ROOT)

print("🔨 Building...")
if subprocess.call(['./install.sh', '-m', 'macos']) != 0:
    print("❌ build failed"); sys.exit(1)

os.environ['EQDYNA2DCYCLEROOT'] = ROOT
os.environ['PATH'] = f"{ROOT}/bin:{ROOT}/scripts:{os.environ['PATH']}"

case = 'work/smoke_case'
subprocess.call(['rm', '-rf', case])
print(f"🏗️  Creating case at {case}...")
if subprocess.call(['python3', 'scripts/create.newcase',
                    '--work_dir', case, '--compset', 'test.subei', '--force']) != 0:
    print("❌ create.newcase failed"); sys.exit(1)

os.chdir(case)
print("🌐 meshgen + case.setup...")
if subprocess.call(['python3', 'meshgen.py']) != 0: sys.exit(1)
if subprocess.call(['python3', 'case.setup']) != 0: sys.exit(1)

print("🚀 run.sh (1 cycle)...")
if subprocess.call(['bash', 'run.sh']) != 0:
    print("❌ run.sh failed"); sys.exit(1)

out = 'aRawSimuData/totalop.txt1'
if not os.path.exists(out) or os.path.getsize(out) == 0:
    print(f"❌ {out} missing or empty"); sys.exit(1)

print(f"✅ smoke test passed ({out} = {os.path.getsize(out)} bytes)")
