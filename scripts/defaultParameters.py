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

    C_mesh = 3
    ntotft = 3
    friclaw = 4 
    dt1 = 0.01 
    term = 200. 
    icstart, icend = 1, 200 
    
    fric_fs = 0.5 
    fric_fd = 0.465
    fric_fv = 0.49
    fric_fini = 0.45
    critd0 = 0.5 # m
    critv0 = 0.2 # m/s 
    critt0 = 0.2 # s
    vrupt0 = 1.5e3 # m/s, forced rupture velocity for nucleateion.
    radius = 2.0e3 # m, radius for forced nucleation patch. 
    
    vp = 6.e3 # m/s
    vs = 3.464e3 # m/s
    rou = 2.67e3 # kg/m^3
    
    eta0 = 8.4e21 # Pa-s, pseudo viscosity for interseismic loading solution.
    # Adjusted according to maximum shear rate, and shear modulus. 
    maxShearStrainLoadRate = 1.427e-14 # On-fault maximum shearing strain loading rate.
    
    ambientnorm = -100.e6 # background ambient on-fault normal stress, Pa. 
    debug = 0 # 1/0, activate/deactivate debugging mode.
    plotmesh = 0 # 1/0, genearte/NOTgenearte mesh files.
    
    yext = 10.e3 # m, external range outside of uniform grid zone along both x & y.
    rat = 1.025 # enlarging ratio for quadralaterals.
    dxy = 300. # m, cell size
    
    ftcn = [15, 10, 80]

    exe = 'run_eqdyna2d_' + _read_version()  # auto-tracks repo VERSION file
    
    