#! /usr/bin/env python3

from defaultParameters import parameters

par = parameters()

par.C_mesh = 3
par.ntotft = 5
par.friclaw = 4 
par.dt1 = 0.01 
par.term = 200. 
par.icstart, par.icend = 1, 10 

par.fric_fs = 0.5 
par.fric_fd = 0.465
par.fric_fv = 0.49
fric_fini = 0.45
par.critd0 = 0.5 # m
par.critv0 = 0.2 # m/s 
par.critt0 = 0.2 # s
par.vrupt0 = 1.5e3 # m/s, forced rupture velocity for nucleateion.
par.radius = 2.0e3 # m, radius for forced nucleation patch. 

par.vp = 6.e3 # m/s
par.vs = 3.464e3 # m/s
par.rou = 2.67e3 # kg/m^3

par.eta0 = 6.e21 # Pa-s, pseudo viscosity for interseismic loading solution.
# Adjusted according to maximum shear rate, and shear modulus. 
par.maxShearStrainLoadRate = 1.427e-14 # On-fault maximum shearing strain loading rate.

par.ambientnorm = -100.e6 # background ambient on-fault normal stress, Pa. 
par.debug = 0 # 1/0, activate/deactivate debugging mode.
par.plotmesh = 0 # 1/0, genearte/NOTgenearte mesh files.

par.yext = 10.e3 # m, external range outside of uniform grid zone along both x & y.
par.rat = 1.025 # enlarging ratio for quadralaterals.
par.dxy = 300. # m, cell size

par.ftcn = [79, 53, 83, 167, 34]
