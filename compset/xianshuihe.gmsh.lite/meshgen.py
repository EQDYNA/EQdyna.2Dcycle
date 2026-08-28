#!/usr/bin/env python3
# meshgen.py — SAF test case
#
# Uses gmsh embedded-curves approach:
#   - One rectangular domain (no sub-surfaces, no connector lines)
#   - Fault lines embedded as interior constraints via gmsh.model.mesh.embed
#   - No shared nodes between non-touching fault traces → no junction slivers
#
# Fault layout (3 independent embedded faults; pre-splined to be smooth):
#   ssaf  (single fault, full SAF trace)
#   sjfn  (independent)
#   sjfs  (independent)

# Xianshuihe fault system, meshed with the EMBEDDED scheme (as in
# saf.gmsh.lite), NOT the chain-of-surfaces scheme in gulang.gmsh.lite.
#
# Every fault is an independent curve embedded in one rectangular surface.
# There are NO connector lines between fault tips. The chain scheme needs
# them, because it uses the faults as the boundary between a top and a bottom
# block, so any gap has to be bridged for the surface loops to close -- but
# nothing in the geology connects one fault tip to the next, and those
# invented lines are what produced the doubling-back hairpin at ft1/ft2, the
# lens surfaces for the splays, and the orphaned split nodes (issue #1).
#
# Embedded faults have free tips, so the overlapping step-overs are carried
# exactly as digitised with no special-casing.
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

ftNamesForGmsh   = ['ft1','ft2','ft3','ft4','ft5','ft6','ft7']
ftNames          = ['ft1','ft2','ft3','ft4','ft5','ft6','ft7']
ftNamesForOutput = ['ft1','ft2','ft3','ft4','ft5','ft6','ft7']

system = "xianshuihe"

dx            = 0.4   # km — target element size on faults (matches par.dxy)
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

# Feed gmsh the original .gmt control points (densifying them broke the
# mesh — gmsh tolerance merges/drops some dense points, leaving orphan
# fault nodes with zero mass). Build a natural cubic spline y(x) from each
# fault's control points; we evaluate its analytical derivative at every
# gmsh fault node later for smooth (C¹) tangents — same intent as
# meshgen1.f90 lines 437-442.
from scipy.interpolate import CubicSpline as _CS
ftCtrlSpline = {}      # {ftName: (xs_km, CubicSpline_obj)}
for ftName in ftNamesForGmsh:
    ftFileName = 'user_fault_geometry_input/' + ftName + '.gmt.txt'
    ftLoc = loadFtLoc(ftFileName)
    ftCtrlSpline[ftName] = (ftLoc[:, 0].copy(), _CS(ftLoc[:, 0], ftLoc[:, 1], bc_type='natural'))
    numOfControlPts, ftEndNodeId, ftRange = createFtPoints(ftLoc, numOfControlPts, dx)
    ftEndNodeIdDict[ftName] = ftEndNodeId
    modelRange = redefineModelRange(modelRange, ftRange)
    lineCount, ftEndLineId = createLinesForFt(ftEndNodeId, lineCount)
    ftEndLineIdDict[ftName] = ftEndLineId

modelRange = extendModelRange(modelRange, ext)

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
ftLinesDict = {k: list(range(ftEndLineIdDict[k][0], ftEndLineIdDict[k][1] + 1))
               for k in ftNamesForGmsh}
allFaultLines = [l for k in ftNamesForGmsh for l in ftLinesDict[k]]
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

# Pure-quad meshing by SUBDIVISION rather than recombination.
#
# gmsh's default path triangulates and then recombines, and where it cannot
# pair triangles it simply leaves them -- fac.txt exports quads only, so those
# elements are dropped, the FE mesh gets holes, and the split nodes bordering
# them orphan (issue #1). Around an embedded fault's free tip the mesher
# closes the front wherever it lands, so that failure is routine there.
#
# SubdivisionAlgorithm 1 subdivides every element into quads instead, so the
# mesh is all-quad BY CONSTRUCTION and recombination cannot fail. Measured on
# this geometry:
#
#   recombination (default)  13511 quads,  8 triangles, 3 orphans, 0 bad, 6 s
#   subdivision              54068 quads,  0 triangles, 0 orphans, 0 bad, 9 s
#
# It also halves the element size, so dx = 0.4 gives ~200 m elements and about
# 4x the cell count. dx = 0.8 with subdivision was tried to recover 400 m and
# is NOT clean -- 28890 quads but 1 badly shaped cell returns.
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.option.setNumber("Mesh.Smoothing", 5)
gmsh.model.mesh.setRecombine(2, 1)
gmsh.model.mesh.generate(2)

os.makedirs('fem_mesh_output', exist_ok=True)
gmsh.write('fem_mesh_output/' + meshName + '.msh')

# ── 6. extract fault nodes from gmsh ────────────────────────────────────────
ftTag = dict(ftLinesDict)

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
# Tangent = analytical derivative of natural cubic spline y(x) built from
# .gmt control points, evaluated at every gmsh fault node's x-coordinate.
# Mirrors meshgen1.f90 lines 437-442 — C¹-smooth across all nodes.
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

# ── 14. fault physics ─────────────────────────────────────────────────────────
ftPhys = defineSysPhys(system, ftNames, xCoorDict, yCoorDict)

# Per-node loading: interpolate from paper.saf.A/Rate_direction.txt
# (proven production loading) onto the gmsh fault nodes by arc-length.
#
# Encoding (matches paper interstress.f90 explicitly):
#   ftLoadMaxShear = γ_per_node                  (= rd(1) per node)
#   ftLoadAngle    = φ_per_node                  (= rd(2) per node)
#   ftLoadWt       = 450  (constant)             (no extra γ encoding via wt)
#   ftVis          = ant0·str / γ_per_node       (= paper's per-node `ant`,
#                                                  inverse-proportional in γ)
# With C_mesh=3 strbld formula `ant = ftVis·450/ftLoadWt`, this yields
# `ant = ftVis = ant0·str/γ`, matching paper's interstress.f90 exactly.
# loading applied post-mesh from GSRM (see strain_rate_loading.py)

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

# Save a paper-Fig.2 style plot of per-node loading inputs for visual sanity check.
plotLoadingInputs(ftPhys, ftNamesForOutput, xCoorDict, yCoorDict)
plotFaults(xCoorDict, yCoorDict, ftPhys, ftNamesForOutput)

if debugMode:
    plotSystemPhys(ftPhys, ftNames, xCoorDict)
