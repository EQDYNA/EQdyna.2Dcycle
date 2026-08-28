#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!pip install matplotlib meshio
import os
import matplotlib
matplotlib.use('Agg')
import gmsh
import numpy as np
import meshio
import matplotlib.pyplot as plt
from meshGenLib import *
from userDefinedFaultSysGeoPhys import *

debugMode = False
meshName = 'eqdynaMesh'

# Xianshuihe through-going chain: 5 faults succeeding one another SE -> NW,
# connected by auxiliary lines for mesh topology. Built by
# export_fault_geometry.py from the 1:1M KML; see README for the mapping
# mf1..mf5 -> ft1,ft2,ft5,ft6,ft7 and which polylines each merges.
#
# The two splay strands (seg3, seg5) are NOT here. They run parallel to mf2
# rather than continuing it, so they need a hand-built multi-surface topology
# like subei.gmsh.lite's, not this chain. Deferred deliberately.
ftNamesForGmsh = ['mf1', 'mf2', 'mf3', 'mf4', 'mf5']
ftNames = ['mf1', 'mf2', 'mf3', 'mf4', 'mf5']


system = "xianshuihe"

# length in km
dx = 0.4
dxAtBoundary = 20
totalSimuTime = 30
vp = 6000
# NOTE: reduced domain for a first look. The reflection-safe extension is
#   ext = 10 + totalSimuTime*vp/1e3   (= 190 km here, 1210 km at term=200 s)
# which keeps boundary reflections from reaching the fault before the rupture
# ends. At ext = 100 km a P-wave round trip is ~33 s, so reflections DO return
# within a 200 s dynamic run. Fine for meshing and geometry checks; restore
# the formula before trusting rupture dynamics.
ext = 100.0

tolerance = 1e-6


# In[2]:


# generating reference mesh with GMSH without split nodes
gmsh.initialize()
#gmsh.model.geo.setFactory("OpenCASCADE")
gmsh.model.add(meshName)
gmsh.option.setNumber("Geometry.Tolerance", 1e-3)
modelRange={'xmin':0,
          'xmax':0,
          'ymin':0,
          'ymax':0}

ftEndNodeIdDict = {}
ftEndLineIdDict = {}
numOfControlPts = 0
lineCount = 0
surfaceCount = 0
from scipy.interpolate import CubicSpline as _CS
ftCtrlSpline = {}      # {ftName: (xs_km, CubicSpline_obj)} for analytical tangent
for ftName in ftNamesForGmsh:
    ftFileName = 'user_fault_geometry_input/' + ftName + '.gmt.txt'
    ftLoc = loadFtLoc(ftFileName)
    # Decimate control points to ~3*mesh-size spacing so gmsh has freedom to
    # place fault nodes uniformly without colliding with control points; the
    # spline derivative then reflects mesh-scale direction (not sub-meter
    # jaggedness from the original .gmt).
    xs_dec, ys_dec = decimateControlPoints(ftLoc[:, 0], ftLoc[:, 1], dx * 3.0)
    ftLoc = np.column_stack([xs_dec, ys_dec])
    ftCtrlSpline[ftName] = (xs_dec.copy(), _CS(xs_dec, ys_dec, bc_type='natural'))
    numOfControlPts, ftEndNodeId, ftRange = createFtPoints(ftLoc, numOfControlPts, dx)
    ftEndNodeIdDict[ftName] = ftEndNodeId
    modelRange = redefineModelRange(modelRange, ftRange)
    lineCount, ftEndLineId = createLinesForFt(ftEndNodeId, lineCount)
    ftEndLineIdDict[ftName] = ftEndLineId
modelRange = extendModelRange(modelRange, ext)


#print(ftRange)
#print(modelRange)
#print(ftEndNodeIdDict)
#print(ftEndLineIdDict)

numOfControlPts, boundaryNodeIdDict = createBoundaryNodes(numOfControlPts, modelRange, dxAtBoundary)
#print(boundaryNodeIdDict)

# linking faults to model boundary
# this part always needs some customized design
# all lines are defined counterclockwisely.

mf1Curve = [i+ftEndLineIdDict['mf1'][0] for i in range(ftEndLineIdDict['mf1'][1]-ftEndLineIdDict['mf1'][0]+1)]
mf2Curve = [i+ftEndLineIdDict['mf2'][0] for i in range(ftEndLineIdDict['mf2'][1]-ftEndLineIdDict['mf2'][0]+1)]
mf3Curve = [i+ftEndLineIdDict['mf3'][0] for i in range(ftEndLineIdDict['mf3'][1]-ftEndLineIdDict['mf3'][0]+1)]
mf4Curve = [i+ftEndLineIdDict['mf4'][0] for i in range(ftEndLineIdDict['mf4'][1]-ftEndLineIdDict['mf4'][0]+1)]
mf5Curve = [i+ftEndLineIdDict['mf5'][0] for i in range(ftEndLineIdDict['mf5'][1]-ftEndLineIdDict['mf5'][0]+1)]
#print(mf1Curve, mf2Curve, mf3Curve, mf4Curve, mf5Curve)

# boundary edges
lineCount += 1
T = gmsh.model.geo.addLine(boundaryNodeIdDict['right_top'], boundaryNodeIdDict['left_top'], tag=lineCount)
lineCount += 1
B = gmsh.model.geo.addLine(boundaryNodeIdDict['right_bottom'], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
L = gmsh.model.geo.addLine(boundaryNodeIdDict['left_top'], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
R = gmsh.model.geo.addLine(boundaryNodeIdDict['right_top'], boundaryNodeIdDict['right_bottom'], tag=lineCount)

# fault chain to boundary
lineCount += 1
FT_LT = gmsh.model.geo.addLine(ftEndNodeIdDict['mf1'][0], boundaryNodeIdDict['left_top'], tag=lineCount) # ft1 start to left top
lineCount += 1
FT_LB = gmsh.model.geo.addLine(ftEndNodeIdDict['mf1'][0], boundaryNodeIdDict['left_bottom'], tag=lineCount) # ft1 start to left bottom
lineCount += 1
FT_RT = gmsh.model.geo.addLine(ftEndNodeIdDict['mf5'][1], boundaryNodeIdDict['right_top'], tag=lineCount) # ft5 end to right top
lineCount += 1
FT_RB = gmsh.model.geo.addLine(ftEndNodeIdDict['mf5'][1], boundaryNodeIdDict['right_bottom'], tag=lineCount) # ft5 end to right bottom

# connecting adjacent faults
lineCount += 1
F1_F2 = gmsh.model.geo.addLine(ftEndNodeIdDict['mf1'][1], ftEndNodeIdDict['mf2'][0], tag=lineCount) # ft1 end to ft2 start
lineCount += 1
F2_F3 = gmsh.model.geo.addLine(ftEndNodeIdDict['mf2'][1], ftEndNodeIdDict['mf3'][0], tag=lineCount) # ft2 end to ft3 start
lineCount += 1
F3_F4 = gmsh.model.geo.addLine(ftEndNodeIdDict['mf3'][1], ftEndNodeIdDict['mf4'][0], tag=lineCount) # ft3 end to ft4 start
lineCount += 1
F4_F5 = gmsh.model.geo.addLine(ftEndNodeIdDict['mf4'][1], ftEndNodeIdDict['mf5'][0], tag=lineCount) # ft4 end to ft5 start

# the fault chain as a single curve list
faultChain = mf1Curve+[F1_F2]+mf2Curve+[F2_F3]+mf3Curve+[F3_F4]+mf4Curve+[F4_F5]+mf5Curve
reversedFaultChain = addMinusToList(mf5Curve[::-1])+[-F4_F5]+addMinusToList(mf4Curve[::-1])+[-F3_F4]+addMinusToList(mf3Curve[::-1])+[-F2_F3]+addMinusToList(mf2Curve[::-1])+[-F1_F2]+addMinusToList(mf1Curve[::-1])

# top block (above fault chain)
surfaceCount = createSurface(faultChain+[FT_RT, T, -FT_LT], surfaceCount)
# left block (triangle)
surfaceCount = createSurface([FT_LB, -L, -FT_LT], surfaceCount)
# bottom block (below fault chain)
surfaceCount = createSurface(reversedFaultChain+[FT_LB, -B, -FT_RB], surfaceCount)
# right block (triangle)
surfaceCount = createSurface([FT_RB, -R, -FT_RT], surfaceCount)

gmsh.model.geo.synchronize()

for iSur in range(surfaceCount):
    surfaceId = iSur+1
    print('recombining surface id ', surfaceId)
    gmsh.model.mesh.setRecombine(2, surfaceId)

#gmsh.model.mesh.coherence() not available in python
gmsh.option.setNumber("Mesh.Smoothing", 2)
gmsh.model.mesh.generate(2)
#gmsh.model.mesh.removeUnusedEntities()
os.makedirs('fem_mesh_output', exist_ok=True)
gmsh.write('fem_mesh_output/' + meshName+'.msh')
#gmsh.finalize()


# In[3]:


ftTag = {}
ftTag['mf1'] = mf1Curve
ftTag['mf2'] = mf2Curve
ftTag['mf3'] = mf3Curve
ftTag['mf4'] = mf4Curve
ftTag['mf5'] = mf5Curve
if debugMode==True:
    print(ftTag)


# In[4]:


nodeTagsFtDict={}
xCoorDict={}
yCoorDict={}
for key in ftNames:
    nodeTagsFtDict[key], xCoorDict[key], yCoorDict[key] = extractFtNodes(ftTag[key])
    if debugMode==True:
        print(nodeTagsFtDict[key])


# In[5]:


# plotting
mesh = meshio.read('fem_mesh_output/' + meshName+'.msh')
points = mesh.points[:, :2]  # Get the x, y coordinates
cells = mesh.cells_dict["quad"]  # Assuming quadrilateral elements

if debugMode==True:
    fig, ax = plt.subplots(figsize=(15, 20), dpi=600)
    ax.scatter(points[:, 0], points[:, 1], s=0.01, color='red', zorder=1)
    for cell in cells:
        vertices = points[cell]
        ax.add_patch(plt.Polygon(vertices, edgecolor='black', linewidth=0.1, fill=False))
    ax.set_aspect('equal')
    plt.savefig('fem_mesh_output/meshWOSplitNode.png', dpi=600)
    plt.close()
    plt.show()


# In[6]:


ftNodeIdsDict={}
for key in ftNames:
    ftNodeIdsDict[key] = locateFtNodeIds(points, xCoorDict[key], yCoorDict[key], tolerance)
    if debugMode==True:
        print(ftNodeIdsDict[key])


# In[7]:


if debugMode==True:
    # check and plot ft nodes
    fig, ax = plt.subplots(figsize=(15, 10), dpi=600)
    ax.scatter(points[:, 0], points[:, 1], s=0.1, color='red', zorder=1)
    for ftName in ftNames:
        ax.scatter(points[ftNodeIdsDict[ftName],0], points[ftNodeIdsDict[ftName],1], s=0.3, color='black', zorder=2)
    plt.savefig('fem_mesh_output/meshWithFaultNodes.png', dpi=600)
    plt.close()
    plt.show()


# In[8]:


elemIdsAboveFtDict = {}
elemIdsBelowFtDict = {}
for key in ftNames:
    elemIdsAboveFtDict[key], elemIdsBelowFtDict[key] = extractIdsforFtElem(ftNodeIdsDict[key], points, cells)

if debugMode==True:
    print(elemIdsAboveFtDict)


# In[9]:


slaveNodeIdsDict={}
masterSlaveNodeIdRelation = {}
pointsWithSplitNodes = np.copy(points)
for ftNameKey in ftNames:
    pointsWithSplitNodes, slaveNodeIdsDict[ftNameKey] = createSplitNodes(ftNodeIdsDict[ftNameKey], pointsWithSplitNodes)

    masterSlaveNodeIdRelation[ftNameKey] = [ftNodeIdsDict[ftNameKey], slaveNodeIdsDict[ftNameKey]]

if debugMode==True:
    print(slaveNodeIdsDict)
    print(ftNodeIdsDict)
    print(' ')
    print(masterSlaveNodeIdRelation)


# In[10]:


if debugMode==True:
    # check and plot ft nodes with split nodes
    fig, ax = plt.subplots(figsize=(15, 10), dpi=600)
    plt.scatter(points[:, 0], points[:, 1], s=0.01, color='red', zorder=1)
    for key in ftNames:
        plt.scatter(pointsWithSplitNodes[ftNodeIdsDict[key],0], pointsWithSplitNodes[ftNodeIdsDict[key],1], s=0.03, color='black', zorder=2)
        plt.scatter(pointsWithSplitNodes[slaveNodeIdsDict[key],0], pointsWithSplitNodes[slaveNodeIdsDict[key],1], marker='*', s=0.003, color='blue', zorder=2)
    plt.close()
    plt.show()


# In[11]:


for ftNameKey in ftNames:
    #print(masterSlaveNodeIdRelation[ftNameKey][1])
    cells = replaceMasterWithSlaveNodes(cells, masterSlaveNodeIdRelation[ftNameKey], elemIdsAboveFtDict[ftNameKey])


# In[12]:


if debugMode==True:
    # testing if node orders are counterclockwise for cell id 1
    vertices = getCellNodeCoors(cells[1], pointsWithSplitNodes)
    showCellNodes(vertices)
    isThisQuadCounterclockwise(vertices)


# In[13]:


## testing if node orders are counterclockwise for cell id 1728
#vertices = getCellNodeCoors(cells[1728], pointsWithSplitNodes)
#showCellNodes(vertices)
#isThisQuadCounterclockwise(vertices)


# In[14]:


reorderedCells = reorderCellNodesCounterclockwise(cells, pointsWithSplitNodes)


# In[15]:


# Tangent = analytical derivative of natural cubic spline y(x) built from
# .gmt control points (matches saf.gmsh.lite + meshgen1.f90). C¹-smooth.
ftNodeTanAndLen = {}
for ftNameKey in ftNames:
    xs_ctrl, cs = ftCtrlSpline[ftNameKey]
    x_node = np.asarray(xCoorDict[ftNameKey], dtype=float)
    y_node = np.asarray(yCoorDict[ftNameKey], dtype=float)
    x_q   = np.clip(x_node, xs_ctrl.min(), xs_ctrl.max())
    slope = cs(x_q, 1)                         # dy/dx, analytical
    mag   = np.sqrt(slope * slope + 1.0)
    tx_n  = 1.0 / mag                          # convention matches meshgen1
    ty_n  = slope / mag
    seg   = np.hypot(np.diff(x_node), np.diff(y_node))   # km
    L     = np.empty(len(x_node))
    L[0]  = 0.5 * seg[0]
    L[-1] = 0.5 * seg[-1]
    L[1:-1] = 0.5 * (seg[:-1] + seg[1:])
    ftNodeTanAndLen[ftNameKey] = [[float(tx_n[i]), float(ty_n[i]), float(L[i])]
                                  for i in range(len(x_node))]


# In[16]:


ftPhys = defineSysPhys(system, ftNames, xCoorDict, yCoorDict)

# Default loading: uniform compression along +x. φ_node = -atan2(ty,tx) (deg).
for ftNameKey in ftNames:
    rows = ftNodeTanAndLen[ftNameKey]
    tx_arr = np.array([r[0] for r in rows])
    ty_arr = np.array([r[1] for r in rows])
    phi    = uniformXLoadingAngle(tx_arr, ty_arr)
    for i in range(len(ftPhys[ftNameKey])):
        ftPhys[ftNameKey][i][3] = float(phi[i])     # ftLoadAngle


# In[17]:


# Change working directory to output folder temporarily
import os
original_dir = os.getcwd()
os.chdir('fem_mesh_output')

meshInfo, pointsWithSplitNodes, reorderedCells, nsmp, nsmpGeoPhys = \
    writeFilesForEQdyna(pointsWithSplitNodes, reorderedCells, masterSlaveNodeIdRelation, ftNodeTanAndLen, ftPhys, modelRange, ftNames)

# Return to original directory
os.chdir(original_dir)


# In[18]:


# Save paper-Fig.2-style plot of per-node loading inputs for sanity check.
plotLoadingInputs(ftPhys, ftNames, xCoorDict, yCoorDict)
plotFaults(xCoorDict, yCoorDict, ftPhys, ftNames)

if debugMode==True:
    plotSystemPhys(ftPhys, ftNames, xCoorDict)
