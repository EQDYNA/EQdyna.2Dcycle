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
# ft1..ft7; see README for which polylines each fault merges.
#
# Seven faults. Five form the through-going chain ft1 -> ft2 -> ft5 -> ft6 ->
# ft7; ft3 and ft4 are the Kangding splay strands. ft3 runs ABOVE the chain
# and ft4 BELOW it, so neither is embedded inside a block -- each divides one
# side into a lens (between the splay and the chain) plus the outer block.
# That gives six surfaces against gulang's four.
ftNamesForGmsh = ['ft1', 'ft2', 'ft3', 'ft4', 'ft5', 'ft6', 'ft7']
ftNames = ['ft1', 'ft2', 'ft3', 'ft4', 'ft5', 'ft6', 'ft7']


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
ftDecimated = {}       # {ftName: decimated control points}
for ftName in ftNamesForGmsh:
    ftFileName = 'user_fault_geometry_input/' + ftName + '.gmt.txt'
    ftLoc = loadFtLoc(ftFileName)
    # Decimate control points to ~3*mesh-size spacing so gmsh has freedom to
    # place fault nodes uniformly without colliding with control points; the
    # spline derivative then reflects mesh-scale direction (not sub-meter
    # jaggedness from the original .gmt).
    # RESAMPLE ALONG THE SPLINE rather than decimating the digitised points.
    #
    # paper.saf.A (C_mesh=2) builds its faults with natural cubic splines in
    # meshgen1.f90, so its nodes lie ON a smooth curve. The gmsh path joins
    # control points with straight lines (createLinesForFt uses addLine), so
    # nodes lie on chords. With the 1:1M digitisation the surviving control
    # points are 2.7-9.1 km apart, far coarser than dx = 400 m, and the
    # polyline then departs from the spline by up to 1185 m on ft2 -- about
    # three element widths.
    #
    # Fitting the spline to the ORIGINAL points and resampling it at ~3*dx
    # keeps the same control-point density gmsh wants while putting those
    # points on the smooth curve, so the chord error drops to the sub-element
    # level and the geometry matches the SAF construction.
    _spline_all = _CS(ftLoc[:, 0], ftLoc[:, 1], bc_type='natural')
    _x0, _x1 = ftLoc[0, 0], ftLoc[-1, 0]
    _n = max(int(np.ceil(abs(_x1 - _x0) / (dx * 3.0))) + 1, 4)
    xs_dec = np.linspace(_x0, _x1, _n)
    ys_dec = _spline_all(xs_dec)
    ftLoc = np.column_stack([xs_dec, ys_dec])
    ftCtrlSpline[ftName] = (xs_dec.copy(), _CS(xs_dec, ys_dec, bc_type='natural'))
    numOfControlPts, ftEndNodeId, ftRange = createFtPoints(ftLoc, numOfControlPts, dx)
    ftEndNodeIdDict[ftName] = ftEndNodeId
    # Keep the decimated control points and the id of the first one. gmsh
    # numbers a fault's points consecutively, so interior point k has node id
    # ftEndNodeId[0] + k -- needed to attach the splay auxiliaries to the
    # nearest point ON the chain instead of to its far end.
    ftDecimated[ftName] = ftLoc
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

ft1Curve = [k+ftEndLineIdDict['ft1'][0] for k in range(ftEndLineIdDict['ft1'][1]-ftEndLineIdDict['ft1'][0]+1)]
ft2Curve = [k+ftEndLineIdDict['ft2'][0] for k in range(ftEndLineIdDict['ft2'][1]-ftEndLineIdDict['ft2'][0]+1)]
ft3Curve = [k+ftEndLineIdDict['ft3'][0] for k in range(ftEndLineIdDict['ft3'][1]-ftEndLineIdDict['ft3'][0]+1)]
ft4Curve = [k+ftEndLineIdDict['ft4'][0] for k in range(ftEndLineIdDict['ft4'][1]-ftEndLineIdDict['ft4'][0]+1)]
ft5Curve = [k+ftEndLineIdDict['ft5'][0] for k in range(ftEndLineIdDict['ft5'][1]-ftEndLineIdDict['ft5'][0]+1)]
ft6Curve = [k+ftEndLineIdDict['ft6'][0] for k in range(ftEndLineIdDict['ft6'][1]-ftEndLineIdDict['ft6'][0]+1)]
ft7Curve = [k+ftEndLineIdDict['ft7'][0] for k in range(ftEndLineIdDict['ft7'][1]-ftEndLineIdDict['ft7'][0]+1)]

# boundary edges
lineCount += 1
T = gmsh.model.geo.addLine(boundaryNodeIdDict['right_top'], boundaryNodeIdDict['left_top'], tag=lineCount)
lineCount += 1
B = gmsh.model.geo.addLine(boundaryNodeIdDict['right_bottom'], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
L = gmsh.model.geo.addLine(boundaryNodeIdDict['left_top'], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
R = gmsh.model.geo.addLine(boundaryNodeIdDict['right_top'], boundaryNodeIdDict['right_bottom'], tag=lineCount)

# chain ends to boundary
lineCount += 1
FT_LT = gmsh.model.geo.addLine(ftEndNodeIdDict['ft1'][0], boundaryNodeIdDict['left_top'], tag=lineCount)
lineCount += 1
FT_LB = gmsh.model.geo.addLine(ftEndNodeIdDict['ft1'][0], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
FT_RT = gmsh.model.geo.addLine(ftEndNodeIdDict['ft7'][1], boundaryNodeIdDict['right_top'], tag=lineCount)
lineCount += 1
FT_RB = gmsh.model.geo.addLine(ftEndNodeIdDict['ft7'][1], boundaryNodeIdDict['right_bottom'], tag=lineCount)

# chain connectors
lineCount += 1
C12 = gmsh.model.geo.addLine(ftEndNodeIdDict['ft1'][1], ftEndNodeIdDict['ft2'][0], tag=lineCount)
lineCount += 1
C25 = gmsh.model.geo.addLine(ftEndNodeIdDict['ft2'][1], ftEndNodeIdDict['ft5'][0], tag=lineCount)
lineCount += 1
C56 = gmsh.model.geo.addLine(ftEndNodeIdDict['ft5'][1], ftEndNodeIdDict['ft6'][0], tag=lineCount)
lineCount += 1
C67 = gmsh.model.geo.addLine(ftEndNodeIdDict['ft6'][1], ftEndNodeIdDict['ft7'][0], tag=lineCount)

# Splay auxiliaries. ft3 hangs above ft2, ft4 below it and on to ft5.
#
# Attach each auxiliary to the control point ON ft2 nearest the splay end,
# not to ft2's far end. A straight line from ft3's end to ft2's end CROSSES
# ft2 (ft2 is still at y ~ -10 under ft3 and only rises to -6.5 at its end),
# which makes the lens loop self-intersecting; gmsh then spins instead of
# failing. Short, near-perpendicular auxiliaries avoid it.
def nearestOnFt2(pt):
    """Index of the ft2 control point nearest pt, and its gmsh node id."""
    k = int(np.argmin(np.hypot(*(ftDecimated['ft2'] - pt).T)))
    return k, ftEndNodeIdDict['ft2'][0] + k

k3s, n3s = nearestOnFt2(ftDecimated['ft3'][0])
k3e, n3e = nearestOnFt2(ftDecimated['ft3'][-1])
k4s, n4s = nearestOnFt2(ftDecimated['ft4'][0])
k3s, k3e = min(k3s, k3e), max(k3s, k3e)
print(f"splay attach points on ft2: ft3 [{k3s}, {k3e}], ft4 start [{k4s}]")

lineCount += 1
A3s = gmsh.model.geo.addLine(ftEndNodeIdDict['ft2'][0] + k3s, ftEndNodeIdDict['ft3'][0], tag=lineCount)
lineCount += 1
A3e = gmsh.model.geo.addLine(ftEndNodeIdDict['ft3'][1], ftEndNodeIdDict['ft2'][0] + k3e, tag=lineCount)
lineCount += 1
A4s = gmsh.model.geo.addLine(ftEndNodeIdDict['ft2'][0] + k4s, ftEndNodeIdDict['ft4'][0], tag=lineCount)
lineCount += 1
A4e = gmsh.model.geo.addLine(ftEndNodeIdDict['ft4'][1], ftEndNodeIdDict['ft5'][0], tag=lineCount)

rev = addMinusToList
# ft2's line ids are consecutive; line j joins control points j and j+1.
ft2Pre  = ft2Curve[:k3s]
ft2Mid  = ft2Curve[k3s:k3e]
ft2Post = ft2Curve[k3e:]
ft2ToF4 = ft2Curve[:k4s]
ft2FromF4 = ft2Curve[k4s:]

# Upper envelope: chain, routed over ft3 for the span the splay covers.
upper = (ft1Curve + [C12] + ft2Pre + [A3s] + ft3Curve + [A3e] + ft2Post
         + [C25] + ft5Curve + [C56] + ft6Curve + [C67] + ft7Curve)
# Lower envelope, ft7 back to ft1, routed under ft4.
lower = (rev(ft7Curve[::-1]) + [-C67] + rev(ft6Curve[::-1]) + [-C56]
         + rev(ft5Curve[::-1]) + [-A4e] + rev(ft4Curve[::-1]) + [-A4s]
         + rev(ft2ToF4[::-1]) + [-C12] + rev(ft1Curve[::-1]))

# top block, above the upper envelope
surfaceCount = createSurface(upper + [FT_RT, T, -FT_LT], surfaceCount)
# lens between ft3 and the ft2 span beneath it
surfaceCount = createSurface([A3s] + ft3Curve + [A3e] + rev(ft2Mid[::-1]), surfaceCount)
# lens between ft4 and the chain above it (ft2 from the attach point, then C25)
surfaceCount = createSurface([A4s] + ft4Curve + [A4e, -C25] + rev(ft2FromF4[::-1]), surfaceCount)
# bottom block, below the lower envelope
surfaceCount = createSurface(lower + [FT_LB, -B, -FT_RB], surfaceCount)
# corner triangles
surfaceCount = createSurface([FT_LB, -L, -FT_LT], surfaceCount)
surfaceCount = createSurface([FT_RB, -R, -FT_RT], surfaceCount)

gmsh.model.geo.synchronize()

for iSur in range(surfaceCount):
    surfaceId = iSur+1
    print('recombining surface id ', surfaceId)
    gmsh.model.mesh.setRecombine(2, surfaceId)

#gmsh.model.mesh.coherence() not available in python

# Mesh size comes from the per-point sizes set in createFtPoints /
# createBoundaryNodes. A gmsh field-based size law is available as
# meshGenLib.applyDistanceSizeField(), OFF here because it measured as a
# regression on this geometry:
#
#   point sizes (this)     10 triangles,  2 orphans, 34602 cells,  5.7 s
#   per-point setSize      14 triangles,  4 orphans, 34866 cells,  ~6 s
#   distance size field    24 triangles,  6 orphans, 57249 cells,  3m27s
#
# All three refinement attempts made recombination WORSE. The only thing that
# cleared every check was trimming the faults apart, reverted because it turns
# the overlapping step-overs into end-to-end gaps. See issue #1.

# NOTE: sizing the overlap control points to the gap width was tried and is a
# no-op here -- the median gap is 0.405 km against dx = 0.400, i.e. already
# 1.01 elements across. Sizing cannot force a single clean row anyway: the
# free mesher places nodes where it likes, and it hugs ft2 (0.07-0.18 km)
# rather than spanning to ft1 (0.30-0.50 km). Forcing one row needs the
# overlap closed into its own 4-sided surface with setTransfiniteSurface.

# Quad-oriented meshing. Without these gmsh uses Algorithm 6 (Frontal-Delaunay
# for TRIANGLES) and RecombinationAlgorithm 1 (simple blossom), which is what
# left stray triangles in the thin ft1/ft2 overlap: fac.txt exports quads only,
# so they were dropped and the bordering split nodes orphaned (issue #1).
#
#   Algorithm 8              Frontal-Delaunay for quads: builds a triangulation
#                            that recombines cleanly
#   RecombinationAlgorithm 3 blossom full-quad: forces full recombination
#   RecombineOptimizeTopology  extra topological clean-up passes
gmsh.option.setNumber("Mesh.Algorithm", 8)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)
# Smoothing 20 / OptimizeTopology 10 was tried and is not better: it just
# trades a >160 deg cell for an aspect>10 one. 2 / 5 is kept.
gmsh.option.setNumber("Mesh.RecombineOptimizeTopology", 5)
gmsh.option.setNumber("Mesh.Smoothing", 2)
gmsh.model.mesh.generate(2)
#gmsh.model.mesh.removeUnusedEntities()
os.makedirs('fem_mesh_output', exist_ok=True)
gmsh.write('fem_mesh_output/' + meshName+'.msh')
#gmsh.finalize()


# In[3]:


ftTag = {}
ftTag['ft1'] = ft1Curve
ftTag['ft2'] = ft2Curve
ftTag['ft3'] = ft3Curve
ftTag['ft4'] = ft4Curve
ftTag['ft5'] = ft5Curve
ftTag['ft6'] = ft6Curve
ftTag['ft7'] = ft7Curve
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
    # Hand a cell back to any split node whose whole fan landed on one side.
    cells, nRepaired = repairOrphanedSplitNodes(
        cells, masterSlaveNodeIdRelation[ftNameKey], pointsWithSplitNodes,
        ftNodeIdsDict[ftNameKey])
    if nRepaired:
        print(f'repaired {nRepaired} orphaned split node(s) on {ftNameKey}')


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
