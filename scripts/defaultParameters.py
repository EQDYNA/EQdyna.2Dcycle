#! /usr/bin/env python3

import os
import numpy as np
from math import *


def _read_version():
    """Read VERSION file at repo root.

    Search order: EQDYNA2DCYCLEROOT env var → walk up from __file__ looking
    for a VERSION file (handles the case where defaultParameters.py was
    copied into a case dir under work/<case>/).
    """
    candidates = []
    root_env = os.environ.get('EQDYNA2DCYCLEROOT')
    if root_env:
        candidates.append(root_env)
    here = os.path.abspath(os.path.dirname(__file__))
    for _ in range(6):
        candidates.append(here)
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    for d in candidates:
        vf = os.path.join(d, 'VERSION')
        if os.path.exists(vf):
            try:
                with open(vf) as f:
                    return f.read().strip()
            except Exception:
                pass
    return 'unknown'


# default from test.tpv1053d
class parameters:

    C_mesh = 3 # mesh mode: 2 = internal Fortran structured mesh from x*_1.txt, 3 = gmsh mesh from meshgen.py
    ntotft = 3 # number of faults
    friclaw = 4 # friction law: 1 slip-weakening, 2 time-weakening, 3 slip-weakening with healing, 4 slip- and rate-weakening (larger of the two), 5 strong rate-weakening
    dt1 = 0.01 # s, time step of the dynamic-rupture solver
    term = 200. # s, maximum duration of one dynamic rupture; an event stops earlier once peak slip rate < 1 mm/s after 5 s
    icstart, icend = 1, 200 # first and last earthquake cycle; icstart > 1 restarts from binaryop
    
    fric_fs = 0.5 # static friction coefficient
    fric_fd = 0.465 # dynamic friction coefficient
    fric_fv = 0.49 # re-strengthening friction coefficient of the rate-weakening branch (friclaw 4)
    fric_fini = 0.45 # initial shear/normal stress ratio on the faults at the start of cycle 1
    critd0 = 0.5 # m, slip-weakening distance Dc
    critv0 = 0.2 # m/s, rate-weakening velocity scale
    critt0 = 0.2 # s, time-weakening duration (friclaw 2 and forced nucleation)
    vrupt0 = 1.5e3 # m/s, rupture velocity imposed inside the nucleation patch
    radius = 2.0e3 # m, radius of the forced nucleation patch
    
    vp = 6.e3 # m/s, P-wave speed
    vs = 3.464e3 # m/s, S-wave speed
    rou = 2.67e3 # kg/m^3, density
    
    eta0 = 8.4e21 # Pa s, reference viscosity of the interseismic loading
    # Adjusted according to maximum shear rate, and shear modulus. 
    maxShearStrainLoadRate = 1.427e-14 # 1/s, reference shear strain rate of the interseismic loading
    
    ambientnorm = -100.e6 # Pa, ambient normal stress on the faults (negative = compression)
    debug = 0 # 1/0, debugging output on/off
    plotmesh = 0 # 1/0, write mesh files for plotting on/off
    # Source-time-function output (stf.bin<icstart>, read with stf_read.py):
    # per event, per slipping fault node, the full dynamic time series of slip
    # rate, slip, shear/normal stress and friction. Read-only: turning it on
    # must not change totalop.txt.
    outputSTF = 0     # 1/0, write per-event rupture time histories (stf.bin) on/off
    stfEvery = 1      # keep every Nth solver step (dt_out = stfEvery*dt1)
    stfVmin = 1.e-3   # m/s; write only nodes whose peak slip rate exceeds this
    
    yext = 10.e3 # m, extent of the coarsening zone outside the uniform mesh, along x and y
    rat = 1.025 # growth ratio of element size in the coarsening zone
    dxy = 300. # m, element size in the uniform mesh zone (C_mesh = 2)
    
    ftcn = [15, 10, 80] # number of geometry control points per fault, in x*_1.txt order (C_mesh = 2)

    exe = 'run_eqdyna2d_' + _read_version()  # solver executable, named from the VERSION file
    
    