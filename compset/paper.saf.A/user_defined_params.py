#! /usr/bin/env python3

from defaultParameters import parameters

par = parameters()

par.C_mesh = 2          # Internal Fortran meshing (structured quad)
par.ntotft = 3
par.friclaw = 4
par.dt1 = 0.01
par.term = 200.0        # paper (Liu et al. 2022) used 200 s for dynamic ruptures
par.icstart, par.icend = 1, 4000

# Paper Model A: fric_fs=0.5, fric_fd=0.465, fric_fv=0.49.
# fric_fini=0.45 is the inverted value that reproduces Pangaea Model A's first interval
# of 540 yr; the archived FE_Global.txt shipped 0.4, which is stale.
par.fric_fs = 0.5
par.fric_fd = 0.465
par.fric_fv = 0.49
par.fric_fini = 0.45
par.critd0 = 0.5        # m
par.critv0 = 0.2        # m/s
par.critt0 = 0.2        # s
par.vrupt0 = 1.5e3      # m/s
par.radius = 2.0e3      # m

par.vp = 6.0e3          # m/s
par.vs = 3.464e3        # m/s
par.rou = 2.67e3        # kg/m^3

par.eta0 = 8.4e21       # Pa-s
par.maxShearStrainLoadRate = 1.427e-14

par.ambientnorm = -100.0e6   # Pa
par.debug = 0
par.plotmesh = 1

par.yext = 10.0e3       # m
par.rat = 1.025
par.dxy = 200.0         # m

par.ftcn = [15, 10, 80]

par.exe = 'run_eqdyna2d_2.0.3'
