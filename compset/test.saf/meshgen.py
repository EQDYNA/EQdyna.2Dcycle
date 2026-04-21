#!/usr/bin/env python3
# meshgen.py — SAF test case
#
# Uses gmsh embedded-curves approach:
#   - One rectangular domain (no sub-surfaces, no connector lines)
#   - Fault lines embedded as interior constraints via gmsh.model.mesh.embed
#   - No shared nodes between non-touching fault traces → no junction slivers
#
# Fault layout:
#   ssaf1 + bridge + ssaf2  →  'ssaf' (one logical fault)
#   sjfn                    →  'sjfn' (independent, not connected to ssaf)
#   sjfs                    →  'sjfs' (independent, not connected to sjfn)

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
meshName  = 'eqdynaMesh'

ftNamesForGmsh = ['ssaf1', 'ssaf2', 'sjfn', 'sjfs']
ftNames        = ['ssaf', 'sjfn', 'sjfs']
ftNamesForOutput = ['sjfn', 'sjfs', 'ssaf']

system = "saf"

dx            = 0.3   # km — target element size on faults
dxAtBoundary  = 20    # km — element size at model boundary
totalSimuTime = 30    # s
vp            = 6000  # m/s
ext           = 10 + totalSimuTime * vp / 1e3   # km — boundary extension

tolerance = 1e-6

# ── 1. gmsh initialisation ──────────────────────────────────────────────────
gmsh.initialize()
gmsh.model.add(meshName)
gmsh.option.setNumber("Geometry.Tolerance", 1e-3)

# ── 2. fault points and lines ───────────────────────────────────────────────
modelRange = {'xmin': 0, 'xmax': 0, 'ymin': 0, 'ymax': 0}
ftEndNodeIdDict = {}
ftEndLineIdDict = {}
numOfControlPts = 0
lineCount       = 0

for ftName in ftNamesForGmsh:
    ftFileName = 'user_fault_geometry_input/' + ftName + '.gmt.txt'
    ftLoc = loadFtLoc(ftFileName)
    numOfControlPts, ftEndNodeId, ftRange = createFtPoints(ftLoc, numOfControlPts, dx)
    ftEndNodeIdDict[ftName] = ftEndNodeId
    modelRange = redefineModelRange(modelRange, ftRange)
    lineCount, ftEndLineId = createLinesForFt(ftEndNodeId, lineCount)
    ftEndLineIdDict[ftName] = ftEndLineId

modelRange = extendModelRange(modelRange, ext)

# Bridge ssaf1-end → ssaf2-start (same physical fault, just split for geometry).
# This is an embedded fault segment — it does NOT connect to sjfn or sjfs.
lineCount += 1
ssafBridgeLine = lineCount
gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf1'][1],
                        ftEndNodeIdDict['ssaf2'][0],
                        tag=ssafBridgeLine)

# ── 3. rectangular domain ───────────────────────────────────────────────────
numOfControlPts += 1; lb = numOfControlPts
gmsh.model.geo.addPoint(modelRange['xmin'], modelRange['ymin'], 0, dxAtBoundary, tag=lb)
numOfControlPts += 1; rb = numOfControlPts
gmsh.model.geo.addPoint(modelRange['xmax'], modelRange['ymin'], 0, dxAtBoundary, tag=rb)
numOfControlPts += 1; rt = numOfControlPts
gmsh.model.geo.addPoint(modelRange['xmax'], modelRange['ymax'], 0, dxAtBoundary, tag=rt)
numOfControlPts += 1; lt = numOfControlPts
gmsh.model.geo.addPoint(modelRange['xmin'], modelRange['ymax'], 0, dxAtBoundary, tag=lt)

lineCount += 1; Bline = lineCount; gmsh.model.geo.addLine(lb, rb, tag=Bline)
lineCount += 1; Rline = lineCount; gmsh.model.geo.addLine(rb, rt, tag=Rline)
lineCount += 1; Tline = lineCount; gmsh.model.geo.addLine(rt, lt, tag=Tline)
lineCount += 1; Lline = lineCount; gmsh.model.geo.addLine(lt, lb, tag=Lline)

gmsh.model.geo.addCurveLoop([Bline, Rline, Tline, Lline], tag=1)
gmsh.model.geo.addPlaneSurface([1], tag=1)

gmsh.model.geo.synchronize()

# ── 4. embed all fault lines into the surface ───────────────────────────────
# Each fault is a free-floating interior constraint; no fault touches another.
ssaf1Lines  = list(range(ftEndLineIdDict['ssaf1'][0], ftEndLineIdDict['ssaf1'][1] + 1))
ssaf2Lines  = list(range(ftEndLineIdDict['ssaf2'][0], ftEndLineIdDict['ssaf2'][1] + 1))
sjfnLines   = list(range(ftEndLineIdDict['sjfn'][0],  ftEndLineIdDict['sjfn'][1]  + 1))
sjfsLines   = list(range(ftEndLineIdDict['sjfs'][0],  ftEndLineIdDict['sjfs'][1]  + 1))

allFaultLines = ssaf1Lines + [ssafBridgeLine] + ssaf2Lines + sjfnLines + sjfsLines
gmsh.model.mesh.embed(1, allFaultLines, 2, 1)

# ── 5. mesh-size field: fine near faults, coarse at boundary ────────────────
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "CurvesList", allFaultLines)
gmsh.model.mesh.field.setNumber(1, "Sampling", 200)

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField",  1)
gmsh.model.mesh.field.setNumber(2, "SizeMin",  dx)
gmsh.model.mesh.field.setNumber(2, "SizeMax",  dxAtBoundary)
gmsh.model.mesh.field.setNumber(2, "DistMin",  dx)
gmsh.model.mesh.field.setNumber(2, "DistMax",  ext * 0.3)

gmsh.model.mesh.field.setAsBackgroundMesh(2)
# Disable other size sources so the field has full control
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.option.setNumber("Mesh.Smoothing", 5)
gmsh.model.mesh.setRecombine(2, 1)
gmsh.model.mesh.generate(2)

os.makedirs('fem_mesh_output', exist_ok=True)
gmsh.write('fem_mesh_output/' + meshName + '.msh')

# ── 6. extract fault nodes from gmsh ────────────────────────────────────────
# ssaf = ssaf1 + bridge + ssaf2 (one continuous logical fault)
ftTag = {}
ftTag['ssaf'] = ssaf1Lines + [ssafBridgeLine] + ssaf2Lines
ftTag['sjfn'] = sjfnLines
ftTag['sjfs'] = sjfsLines

nodeTagsFtDict = {}
xCoorDict      = {}
yCoorDict      = {}
for key in ftNames:
    nodeTagsFtDict[key], xCoorDict[key], yCoorDict[key] = extractFtNodes(ftTag[key])

# ── 7. read mesh with meshio ─────────────────────────────────────────────────
mesh   = meshio.read('fem_mesh_output/' + meshName + '.msh')
points = mesh.points[:, :2]
cells  = mesh.cells_dict["quad"]

if debugMode:
    fig, ax = plt.subplots(figsize=(15, 20), dpi=300)
    ax.scatter(points[:, 0], points[:, 1], s=0.01, color='red', zorder=1)
    for cell in cells:
        vertices = points[cell]
        ax.add_patch(plt.Polygon(vertices, edgecolor='black', linewidth=0.1, fill=False))
    ax.set_aspect('equal')
    plt.savefig('fem_mesh_output/meshWOSplitNode.png', dpi=300)
    plt.close()

# ── 8. locate fault nodes in meshio node list ────────────────────────────────
ftNodeIdsDict = {}
for key in ftNames:
    print(f'locateFtNodeIds: {key} ...', flush=True)
    ftNodeIdsDict[key] = locateFtNodeIds(points, xCoorDict[key], yCoorDict[key], tolerance)
    print(f'locateFtNodeIds: {key} done', flush=True)

# ── 9. classify elements as above / below each fault ─────────────────────────
elemIdsAboveFtDict = {}
elemIdsBelowFtDict = {}
for key in ftNames:
    elemIdsAboveFtDict[key], elemIdsBelowFtDict[key] = \
        extractIdsforFtElem(ftNodeIdsDict[key], points, cells)

if debugMode:
    fig, ax = plt.subplots(figsize=(15, 10), dpi=300)
    ax.scatter(points[:, 0], points[:, 1], s=0.1, color='red', zorder=1)
    for ftName in ftNames:
        ax.scatter(points[ftNodeIdsDict[ftName], 0],
                   points[ftNodeIdsDict[ftName], 1], s=0.3, color='black', zorder=2)
    plt.savefig('fem_mesh_output/meshWithFaultNodes.png', dpi=300)
    plt.close()

# ── 10. create split nodes ────────────────────────────────────────────────────
slaveNodeIdsDict         = {}
masterSlaveNodeIdRelation = {}
pointsWithSplitNodes     = np.copy(points)
for ftNameKey in ftNames:
    pointsWithSplitNodes, slaveNodeIdsDict[ftNameKey] = \
        createSplitNodes(ftNodeIdsDict[ftNameKey], pointsWithSplitNodes)
    masterSlaveNodeIdRelation[ftNameKey] = \
        [ftNodeIdsDict[ftNameKey], slaveNodeIdsDict[ftNameKey]]

# ── 11. replace master with slave nodes in above-fault elements ───────────────
for ftNameKey in ftNames:
    cells = replaceMasterWithSlaveNodes(
        cells, masterSlaveNodeIdRelation[ftNameKey], elemIdsAboveFtDict[ftNameKey])

# ── 12. reorder cell nodes counter-clockwise ──────────────────────────────────
reorderedCells = reorderCellNodesCounterclockwise(cells, pointsWithSplitNodes)

# ── 13. fault tangent and length ──────────────────────────────────────────────
ftNodeTanAndLen = {}
for ftNameKey in ftNames:
    ftNodeTanAndLen[ftNameKey] = calcTanAndLen(xCoorDict[ftNameKey], yCoorDict[ftNameKey])

# ── 14. fault physics ─────────────────────────────────────────────────────────
ftPhys = defineSysPhys(system, ftNames, xCoorDict, yCoorDict)

# Per-node loading interpolation from 4-column gmt files
REFERENCE_LOAD_RATE = 1.427e-14  # s^-1, equivalent to the legacy 450 nrad/yr reference
REFERENCE_LOAD_WEIGHT = 450.0

def _arc_lengths(xy):
    d = np.sqrt(np.sum(np.diff(xy, axis=0)**2, axis=1))
    return np.concatenate([[0.0], np.cumsum(d)])

_gmt_map = {
    'ssaf': ['user_fault_geometry_input/ssaf1.gmt.txt',
             'user_fault_geometry_input/ssaf2.gmt.txt'],
    'sjfn': ['user_fault_geometry_input/sjfn.gmt.txt'],
    'sjfs': ['user_fault_geometry_input/sjfs.gmt.txt'],
}
for ftNameKey in ftNames:
    ctrl_parts   = [np.loadtxt(f) for f in _gmt_map[ftNameKey]]
    ctrl         = np.vstack(ctrl_parts)          # (N,4): x(km), y(km), rate, angle
    ctrl_s       = _arc_lengths(ctrl[:, :2])
    ft_xy        = np.column_stack([xCoorDict[ftNameKey], yCoorDict[ftNameKey]])
    node_s       = _arc_lengths(ft_xy)
    rate_interp  = np.interp(node_s, ctrl_s, ctrl[:, 2])
    angle_interp = np.interp(node_s, ctrl_s, ctrl[:, 3])
    weight_interp = REFERENCE_LOAD_WEIGHT * rate_interp / REFERENCE_LOAD_RATE
    for i in range(len(ftPhys[ftNameKey])):
        ftPhys[ftNameKey][i][4] = float(weight_interp[i])
        ftPhys[ftNameKey][i][3] = float(angle_interp[i])

# ── 15. write EQdyna input files ──────────────────────────────────────────────
original_dir = os.getcwd()
os.chdir('fem_mesh_output')

meshInfo, pointsWithSplitNodes, reorderedCells, nsmp, nsmpGeoPhys = \
    writeFilesForEQdyna(
        pointsWithSplitNodes,
        reorderedCells,
        masterSlaveNodeIdRelation,
        ftNodeTanAndLen,
        ftPhys,
        modelRange,
        ftNamesForOutput,
    )

os.chdir(original_dir)

if debugMode:
    plotSystemPhys(ftPhys, ftNames, xCoorDict)
