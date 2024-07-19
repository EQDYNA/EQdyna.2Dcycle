#! /usr/bin/env python3

import numpy as np
from math import *

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
    
    eta0 = 6.e21 # Pa-s, pseudo viscosity for interseismic loading solution.
    # Adjusted according to maximum shear rate, and shear modulus. 
    maxShearStrainLoadRate = 1.427e-14 # On-fault maximum shearing strain loading rate.
    
    ambientnorm = -100.e6 # background ambient on-fault normal stress, Pa. 
    debug = 0 # 1/0, activate/deactivate debugging mode.
    plotmesh = 0 # 1/0, genearte/NOTgenearte mesh files.
    
    yext = 10.e3 # m, external range outside of uniform grid zone along both x & y.
    rat = 1.025 # enlarging ratio for quadralaterals.
    dxy = 300. # m, cell size
    
    ftcn = [15, 10, 80]
    
    