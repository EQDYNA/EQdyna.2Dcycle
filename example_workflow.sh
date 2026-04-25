#!/usr/bin/env bash
# example_workflow.sh — minimal end-to-end run of compset paper.saf.A
#
# Reproduces Liu et al. (2022) JGR Solid Earth Model A (Southern San Andreas).
# C_mesh=2 (Fortran structured quad), 4000 cycles, ~200 s dynamic-rupture window
# per cycle. Default is 1 OMP thread (set OMP_NUM_THREADS=8 or 16 before
# running for faster wall-clock); kill with Ctrl-C
# whenever you have enough cycles for what you're testing.
#
# Prereqs: gfortran, python3 (numpy, matplotlib, xarray).
# Run from the repo root.

set -euo pipefail

# 1. One-time install: builds bin/run_eqdyna2d_<VERSION> and exports
#    EQDYNA2DCYCLEROOT + PATH for bin/ and scripts/. Sourced so env
#    survives into the create.newcase / case.setup steps below.
source ./install.sh -m ubuntu          # use -m macos on a Mac, or -e ubuntu / -e macos to also install python deps

# 2. Create a fresh work dir from the compset.
#
#    To use a different compset, change CASE and COMPSET below. List
#    available compsets with:  create.newcase --list
#    Common choices:
#      COMPSET=paper.saf.A         # C_mesh=2 (Fortran mesh), Liu et al. 2022 Model A — DEFAULT
#      COMPSET=saf.gmsh.lite       # C_mesh=3 (gmsh mesh), paper-A reproduction
#      COMPSET=subei.gmsh.lite     # C_mesh=3, Subei fault system (atf, dxs, sbt)
#      COMPSET=gulang.gmsh.lite    # C_mesh=3, Gulang 5-fault system
#    For shorter runs, edit user_defined_params.py inside $CASE after
#    create.newcase to lower icend (default 4000).
CASE=work/paper.saf.A.demo
COMPSET=paper.saf.A
create.newcase --work_dir "$CASE" --compset "$COMPSET" --force
cd "$CASE"

# 3. C_mesh=3 compsets (the .gmsh.lite ones) need an extra mesh step here:
#      python3 meshgen.py
#    paper.saf.A is C_mesh=2 so its mesh is built inside Fortran — skip.

# 4. Generate FE_*.txt + run.sh from user_defined_params.py.
python3 case.setup

# 5. Launch. run.sh nohups the binary in the background, logs to
#    run_<timestamp>.log, and moves outputs into aRawSimuData/ when done.
#    OMP_NUM_THREADS defaults to 1 (serial); override before this line if needed.
bash run.sh

echo
echo "Logs:    $CASE/run_*.log"
echo "Outputs: $CASE/aRawSimuData/{cyclelog,interval,totalop}.txt*"
echo "Plots:   plotRuptureDynamics  (run from inside $CASE after some cycles)"
