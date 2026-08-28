#! /usr/env/bin python3
import numpy as np
import gmsh
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def plotLoadingInputs(ftPhys, ftNames, xCoorDict, yCoorDict, out_path='aPlots/loading_inputs.png'):
    """Plot per-fault loading inputs (γ, φ, ftLoadWt, ftVis) along arc-length.

    Paper Fig. 2 style. Called at end of meshgen.py to verify the per-node
    loading written to nsmpGeoPhys.txt looks right (smooth, paper-like).
    """
    import os
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    colors = ['#1f77b4', '#d62728', '#222222', '#2ca02c']
    fig, ax = plt.subplots(4, 1, figsize=(10, 12), sharex=False)

    for ift, ftName in enumerate(ftNames):
        rows = np.asarray(ftPhys[ftName])  # (n, 6) = [type, dip, γ, φ, wt, vis]
        x_km = np.asarray(xCoorDict[ftName])
        y_km = np.asarray(yCoorDict[ftName])
        seg = np.hypot(np.diff(x_km), np.diff(y_km))
        s_km = np.concatenate(([0.0], np.cumsum(seg)))

        c = colors[ift % len(colors)]
        ax[0].plot(s_km, rows[:, 2], color=c, label=f'{ftName} (n={len(rows)})')
        ax[1].plot(s_km, rows[:, 3], color=c)
        ax[2].plot(s_km, rows[:, 4], color=c)
        ax[3].plot(s_km, rows[:, 5], color=c)

    ax[0].set_ylabel('ftLoadMaxShear (s$^{-1}$)')
    ax[1].set_ylabel('ftLoadAngle φ (deg)')
    ax[2].set_ylabel('ftLoadWt (= 450·γ/str)')
    ax[3].set_ylabel('ftVis (Pa·s)')
    ax[3].set_xlabel('arc-length along fault (km)')
    ax[0].legend(loc='best', fontsize=8)
    ax[0].set_title('Per-node loading inputs written to nsmpGeoPhys.txt')

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def splineFaultFromControl(xs_km, ys_km, dxy_m):
    """Mirror meshgen1.f90: natural cubic spline y(x) through .gmt control
    points, sampled at uniform x spacing dxy_m.  Returns
    (x_dense_km, y_dense_km, tx, ty) with y on the smooth curve and
    tx, ty = analytical spline derivative at each sample.

    Requires xs to be strictly monotonic (same constraint as meshgen1).
    """
    xs = np.asarray(xs_km, dtype=float)
    ys = np.asarray(ys_km, dtype=float)
    cs = CubicSpline(xs, ys, bc_type='natural')
    n = max(2, int(round((xs[-1] - xs[0]) * 1.0e3 / dxy_m)) + 1)
    x_dense = np.linspace(xs[0], xs[-1], n)
    y_dense = cs(x_dense)                  # y on the smooth curve
    slope   = cs(x_dense, 1)               # dy/dx, analytical
    mag     = np.sqrt(slope * slope + 1.0)
    tx      = 1.0 / mag                    # convention matches meshgen1 (line 441)
    ty      = slope / mag
    return x_dense, y_dense, tx, ty


def plotFaults(xCoorDict, yCoorDict, ftPhys, ftNames, out_path='aPlots/faults.png',
               loading_arrow=True):
    """Plot fault traces in map view (x-y, km), colour-coded by ftType, with
    optional loading-direction arrow (+x = along x).

    Saves a PNG showing fault geometry as a sanity check before running.
    """
    import os
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    type_colors = {1: '#1f77b4', -1: '#d62728', 2: '#ff7f0e', -2: '#2ca02c'}
    type_labels = {1: 'left-strike (ftType=1)', -1: 'right-strike (ftType=-1)',
                   2: 'thrust (ftType=2)', -2: 'normal (ftType=-2)'}

    fig, ax = plt.subplots(figsize=(11, 7))
    seen_types = set()
    for ftName in ftNames:
        x = np.asarray(xCoorDict[ftName], dtype=float)
        y = np.asarray(yCoorDict[ftName], dtype=float)
        ftType = int(ftPhys[ftName][0][0])
        c = type_colors.get(ftType, 'k')
        label = type_labels.get(ftType, f'ftType={ftType}') if ftType not in seen_types else None
        seen_types.add(ftType)
        ax.plot(x, y, '-', color=c, linewidth=2.0, label=label)
        # endpoint markers + name
        ax.plot([x[0], x[-1]], [y[0], y[-1]], 'o', color=c, markersize=4)
        ax.annotate(ftName, (x[len(x)//2], y[len(y)//2]),
                    fontsize=8, color=c, ha='center', va='bottom',
                    xytext=(0, 4), textcoords='offset points')

    if loading_arrow:
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        arrow_y = y_min + 0.05 * (y_max - y_min)
        arrow_len = 0.15 * (x_max - x_min)
        ax.annotate('', xy=(x_min + 0.30 * (x_max - x_min) + arrow_len, arrow_y),
                    xytext=(x_min + 0.30 * (x_max - x_min), arrow_y),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=2))
        ax.text(x_min + 0.30 * (x_max - x_min) + arrow_len / 2, arrow_y - 0.04 * (y_max - y_min),
                'loading: max-shear along +x', color='gray', ha='center', fontsize=9)

    ax.set_aspect('equal')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('Fault layout')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def decimateControlPoints(xs_km, ys_km, dx_km):
    """Keep first + last control points; thin intermediate points so consecutive
    kept points are at least dx_km apart in arc-length. Avoids over-resolving
    fault traces that ship with dense (e.g. 1 m) control points, which would
    otherwise let the spline derivative track sub-mesh-scale jaggedness."""
    xs = np.asarray(xs_km, dtype=float)
    ys = np.asarray(ys_km, dtype=float)
    keep = [0]
    last_x, last_y = xs[0], ys[0]
    for i in range(1, len(xs) - 1):
        if np.hypot(xs[i] - last_x, ys[i] - last_y) >= dx_km:
            keep.append(i)
            last_x, last_y = xs[i], ys[i]
    keep.append(len(xs) - 1)
    return xs[keep], ys[keep]


def uniformXLoadingAngle(tx, ty):
    """Default loading model: uniform max-shear strain rate along +x axis.

    Paper convention: φ = fault_strike − max_shear_direction (degrees).
    With max-shear along +x (= 0°), φ = α = atan2(ty, tx) per node.
    Then `rs = γ·cos(2φ)·ant`:
        horizontal fault (α=0):  cos(0)   = 1   → max positive shear
        45° fault       (α=45):  cos(90°) = 0   → no shear
        vertical fault  (α=90):  cos(180°) = -1 → max negative shear

    The existing upper-clamp at 45° in interstress.f90 protects vertical
    faults from runaway negative loading by capping any |α|>45° to 45°.
    """
    return np.degrees(np.arctan2(np.asarray(ty, dtype=float),
                                  np.asarray(tx, dtype=float)))


def loadPaperRateDirAlongArcLen(paper_root, fault_idx, n_active):
    """Reconstruct paper.saf.A's per-fault arc-length and (γ, φ) by re-running
    the meshgen1.f90-style natural cubic spline through control points
    x{fault_idx}_1.txt, then sampling at n_active uniform x-positions.

    Returns (s_paper [m], gamma [rad/s], phi [deg]) all length n_active.
    Used to interpolate paper Rate_direction.txt loading onto a gmsh fault.

    fault_idx: 1-based (1=fault 1, 2=fault 2, 3=fault 3 in paper indexing).
    n_active: number of active rows for this fault in Rate_direction.txt.
    """
    import os
    ctrl = np.loadtxt(os.path.join(paper_root, f'x{fault_idx}_1.txt'))
    xs_km, ys_km = ctrl[:, 0], ctrl[:, 1]
    xs, ys = xs_km * 1.0e3, ys_km * 1.0e3  # m, matches meshgen1

    cs = CubicSpline(xs, ys, bc_type='natural')
    x_paper = np.linspace(xs.min(), xs.max(), n_active)
    y_paper = cs(x_paper)
    seg = np.hypot(np.diff(x_paper), np.diff(y_paper))
    s_paper = np.concatenate([[0.0], np.cumsum(seg)])

    rd = np.loadtxt(os.path.join(paper_root, 'Rate_direction.txt'))
    n_rows = len(rd)
    ntotft = 3
    maxftnode = n_rows // ntotft
    start = (fault_idx - 1) * maxftnode
    rd_fault = rd[start:start + n_active]
    gamma = rd_fault[:, 0]
    phi   = rd_fault[:, 1]
    return s_paper, gamma, phi


def loadFtLoc(ftName):
    
    with open(ftName, 'r') as f:
        lines = f.readlines()
        ftLoc = []
        for line in lines[:]:
            loc = [float(var) for var in line.split()]
            ftLoc.append(loc)

        ftLocUTM = np.array(ftLoc)
    return ftLocUTM

def createFtPoints(ftLocUTM, numOfControlPts, dx):
    ftEndNodeId=np.zeros(2, dtype=int)
    ftRange = {'xmin':ftLocUTM[:,0].min(), 
              'xmax':ftLocUTM[:,0].max(), 
              'ymin':ftLocUTM[:,1].min(), 
              'ymax':ftLocUTM[:,1].max()}
    
    for i in range(ftLocUTM.shape[0]):
        numOfControlPts += 1
        gmsh.model.geo.addPoint(ftLocUTM[i, 0], ftLocUTM[i, 1], 0, meshSize=dx, tag=numOfControlPts)
        if i == 0:
            ftEndNodeId[0] = numOfControlPts
    
    ftEndNodeId[1] = numOfControlPts
    print(numOfControlPts,' nodes have been created so far.')
    print('Fault ends node Ids are ', ftEndNodeId)
    return numOfControlPts, ftEndNodeId, ftRange

def redefineModelRange(modelRange, ftRange):
    if modelRange['xmin']>ftRange['xmin']:
        modelRange['xmin'] = ftRange['xmin']

    if modelRange['ymin']>ftRange['ymin']:
        modelRange['ymin'] = ftRange['ymin']
        
    if modelRange['xmax']<ftRange['xmax']:
        modelRange['xmax'] = ftRange['xmax']

    if modelRange['ymax']<ftRange['ymax']:
        modelRange['ymax'] = ftRange['ymax']
    
    return modelRange

def extendModelRange(modelRange, ext):
    modelRange['xmin'] = modelRange['xmin'] - ext
    modelRange['ymin'] = modelRange['ymin'] - ext
    modelRange['xmax'] = modelRange['xmax'] + ext
    modelRange['ymax'] = modelRange['ymax'] + ext
    return modelRange

def createBoundaryNodes(numOfControlPts, modelRange, dxAtBoundary):
    boundaryNodeIdDict = {} #np.zeros(4, dtype=int)
    
    numOfControlPts+=1
    gmsh.model.geo.addPoint(modelRange['xmin'], modelRange['ymin'], 0, meshSize=dxAtBoundary, tag=numOfControlPts)
    boundaryNodeIdDict['left_bottom'] = numOfControlPts
    
    numOfControlPts+=1
    gmsh.model.geo.addPoint(modelRange['xmax'], modelRange['ymin'], 0, meshSize=dxAtBoundary, tag=numOfControlPts)
    boundaryNodeIdDict['right_bottom'] = numOfControlPts
    
    numOfControlPts+=1
    gmsh.model.geo.addPoint(modelRange['xmax'], modelRange['ymax'], 0, meshSize=dxAtBoundary, tag=numOfControlPts)
    boundaryNodeIdDict['right_top'] = numOfControlPts
    
    numOfControlPts+=1
    gmsh.model.geo.addPoint(modelRange['xmin'], modelRange['ymax'], 0, meshSize=dxAtBoundary, tag=numOfControlPts)
    boundaryNodeIdDict['left_top'] = numOfControlPts
    
    return numOfControlPts, boundaryNodeIdDict

def createLinesForFt(ftEndNodeId, lineCount):
    ftEndLineId = np.zeros(2, dtype=int)
    
    for i in range(ftEndNodeId[1]-ftEndNodeId[0]):
        nodeId = i+ftEndNodeId[0]
        lineCount += 1 
        gmsh.model.geo.addLine(nodeId, nodeId+1, tag=lineCount)
        if i == 0:
            ftEndLineId[0] = lineCount
    
    ftEndLineId[1] = lineCount
    print(lineCount, ' lines have been created so far.')
    print('Ft ends line ids are', ftEndLineId)
    return lineCount, ftEndLineId

def createSurface(curveLoop, surfaceCount):
    surfaceCount += 1
    gmsh.model.geo.addCurveLoop(curveLoop, tag=surfaceCount)
    domain = gmsh.model.geo.addPlaneSurface([surfaceCount], tag=surfaceCount)
    return surfaceCount

def addMinusToList(aList):
    return [-x for x in aList]

def extractFtNodes(ftTag):
    nodeTagsFt = []
    nodeCoorsFt = []
    for lineTag in ftTag:
        nodeTags, nodeCoors, _ = gmsh.model.mesh.getNodes(1, lineTag, True)
        nodeTagsFt += nodeTags.tolist()
        nodeCoorsFt += nodeCoors.tolist()
    
    xCoors = np.array(nodeCoorsFt[0::3])
    yCoors = np.array(nodeCoorsFt[1::3])
    nodeTagsFtArr = np.array(nodeTagsFt)
    
    uniquenodeTagsFt, uniqueIndices = np.unique(nodeTagsFtArr, return_index=True)
    uniqueXCoors = xCoors[uniqueIndices]
    uniqueYCoors = yCoors[uniqueIndices]
    
    sortedIndices = np.argsort(uniqueXCoors)
    sortedXCoors = uniqueXCoors[sortedIndices]
    sortedYCoors = uniqueYCoors[sortedIndices]
    sortednodeTagsFt = uniquenodeTagsFt[sortedIndices]
    #sortedIndices = sorted(range(len(xCoors)), key=lambda i: xCoors[i])
    #sortedXCoors = [xCoors[i] for i in sortedIndices]
    #sortedYCoors = [yCoors[i] for i in sortedIndices]
    #sortednodeTagsFt = [nodeTagsFt[i] for i in sortedIndices]
    #print(len(xCoors), len(nodeTagsFt))
    
    #return nodeTagsFt, xCoors, yCoors
    return list(sortednodeTagsFt), list(sortedXCoors), list(sortedYCoors)

def calcCenterLoc(coors):
    xCoors = [pt[0] for pt in coors]
    yCoors = [pt[1] for pt in coors]
    
    xCenterCoor = sum(xCoors)/len(xCoors)
    yCenterCoor = sum(yCoors)/len(yCoors)
    return [xCenterCoor, yCenterCoor]

def judgeElemDirect(cxOnFt, cxOffFt, ftTangent=None):
    """Which side of the fault a cell lies on. 1/2 = above, 3/4 = below.

    With ftTangent, the side is the sign of the cross product of the local
    fault tangent with the vector to the off-fault centroid -- the actual
    fault normal. Fault points are ordered along +x, so "left of tangent"
    is the +y side and the result agrees with the old global-quadrant test
    wherever the fault is gently sloping.

    The global test (sign of dy alone) is WRONG for a steeply running
    fault: cells on BOTH sides can then have dy > 0, so both are called
    "above", the master node loses its whole cell fan, and it orphans.
    That produced exactly one orphaned master on the steep NE end of the
    Xianshuihe ft4. ftTangent=None keeps the old behaviour for callers
    that have not been updated.
    """
    elemDirectVec = np.array(cxOffFt) - np.array(cxOnFt)
    if ftTangent is not None:
        t = np.asarray(ftTangent, dtype=float)
        nrm = np.hypot(t[0], t[1])
        if nrm > 0:
            cross = t[0]*elemDirectVec[1] - t[1]*elemDirectVec[0]
            if cross > 0:
                return 1          # left of tangent -> above
            elif cross < 0:
                return 3          # right of tangent -> below
            # cross == 0 (degenerate) falls through to the quadrant test
    quadrant = 0
    if elemDirectVec[0]>0 and elemDirectVec[1]>0:
        quadrant = 1
    elif elemDirectVec[0]<0 and elemDirectVec[1]>0:
        quadrant = 2
    elif elemDirectVec[0]<0 and elemDirectVec[1]<0:
        quadrant = 3
    elif elemDirectVec[0]>0 and elemDirectVec[1]<0:
        quadrant = 4
    
    if quadrant==0: 
        print('Error judgement of quandrant and the elem vector is problematic.')
    
    return quadrant

# Function to compare two points with a tolerance
def isTwoPtsClose(point1, point2, tol):
    return all(abs(a - b) < tol for a, b in zip(point1, point2))

def locateFtNodeIds(points, xCoors, yCoors, tolerance):
    pointsArr = np.array(points)
    ftNodeIds = []
    for xcoor, ycoor in zip(xCoors, yCoors):
        target = np.array([xcoor, ycoor])
        for index, point in enumerate(points):
            if isTwoPtsClose(point, target, tolerance):
                #print(f"The index of the target coordinates is: {index}")
                break
        else:
            print("The target coordinates are not in the list.")
        #print(target, ' is found, and node index is ', index) 
        ftNodeIds += [index]
    return ftNodeIds

#def isThisNodeInCell(nodeId, nodeIdsInCell):
#    if 
def extractIdsforFtElem(ftNodeIds, points, cells):
    # Build node->elements lookup once: O(nElems) instead of O(nFtNodes*nElems)
    nodeToElems = {}
    for iElem, nodeIds in enumerate(cells):
        for nid in nodeIds:
            if nid not in nodeToElems:
                nodeToElems[nid] = []
            nodeToElems[nid].append(iElem)

    elemIdsAboveFt = []
    elemIdsBelowFt = []
    classified = set()

    # Pass 1: edge-on-fault cells (2 consecutive fault nodes in the cell)
    for iNode in range(len(ftNodeIds)-1):
        n0, n1 = ftNodeIds[iNode], ftNodeIds[iNode+1]
        candidates = set(nodeToElems.get(n0, [])) & set(nodeToElems.get(n1, []))
        for iElem in candidates:
            if iElem in classified:
                continue
            nodeIds = cells[iElem]
            otherTwoNodes = [tag for tag in nodeIds if tag not in (n0, n1)]
            coorsOnFt  = [points[n0, :], points[n1, :]]
            coorsOffFt = [points[tag, :] for tag in otherTwoNodes]
            cxOnFt  = calcCenterLoc(coorsOnFt)
            cxOffFt = calcCenterLoc(coorsOffFt)
            quadrant = judgeElemDirect(cxOnFt, cxOffFt,
                                       ftTangent=points[n1, :] - points[n0, :])
            if quadrant == 1 or quadrant == 2:
                elemIdsAboveFt.append(iElem)
            if quadrant == 3 or quadrant == 4:
                elemIdsBelowFt.append(iElem)
            classified.add(iElem)

    # Pass 2: corner-only cells (exactly one fault node in the cell, no edge on fault)
    ftSet = set(ftNodeIds)
    for iFt in ftNodeIds:
        for iElem in nodeToElems.get(iFt, []):
            if iElem in classified:
                continue
            nodeIds = cells[iElem]
            ftInCell = [nid for nid in nodeIds if nid in ftSet]
            if len(ftInCell) != 1:
                continue  # skip cells with 0 or ≥2 fault nodes (handled in pass 1 or not relevant)
            otherThree = [tag for tag in nodeIds if tag != iFt]
            cxOnFt  = list(points[iFt, :])
            cxOffFt = calcCenterLoc([points[tag, :] for tag in otherThree])
            # local tangent from this fault node's neighbours along the fault
            kFt = ftNodeIds.index(iFt) if isinstance(ftNodeIds, list) else int(np.where(np.asarray(ftNodeIds) == iFt)[0][0])
            kPrev = max(kFt - 1, 0)
            kNext = min(kFt + 1, len(ftNodeIds) - 1)
            tang = points[ftNodeIds[kNext], :] - points[ftNodeIds[kPrev], :]
            quadrant = judgeElemDirect(cxOnFt, cxOffFt, ftTangent=tang)
            if quadrant == 1 or quadrant == 2:
                elemIdsAboveFt.append(iElem)
            if quadrant == 3 or quadrant == 4:
                elemIdsBelowFt.append(iElem)
            classified.add(iElem)

    return elemIdsAboveFt, elemIdsBelowFt

def createSplitNodes(ftNodeIds, points):
    countSlaveNode = 0
    slaveNodeIds = []
    slaveNodeCoors = []
    totalNodes = points.shape[0]
    for nodeId in ftNodeIds:
        masterNodeCoors = points[nodeId,:]
        #print('coors by getNode', masterNodeCoors)
        slaveNodeCoors = points[nodeId,:]
        points = np.vstack((points, slaveNodeCoors))
        slaveNodeIds += [len(points)-1]
    return points, slaveNodeIds

def replaceMasterWithSlaveNodes(cells, masterSlaveNodeIdRelation, elemIdsAboveFt):
    for i, nodeIdsInCell in enumerate(cells[elemIdsAboveFt]):
        #print('Processing cell, before replacement ', cell)
        for j in range(4):
            try:
                index = masterSlaveNodeIdRelation[0].index(nodeIdsInCell[j])
                cells[elemIdsAboveFt[i]][j] = masterSlaveNodeIdRelation[1][index]
            except:
                index = None
                #print('Skipping this node ', cell[i])
        #print('After replacement ', cell)

    return cells


def getCellNodeCoors(nodeIds, points):
    vertices = np.zeros((4,2))
    for i in range(4):
        vertices[i,:] = points[nodeIds[i],:]
    return vertices

def showCellNodes(vertices):
    fig, ax = plt.subplots()
    label = ['1','2','3','4']
    for i in range(4):
        ax.scatter(vertices[i,0], vertices[i,1], label=label[i])
    
    ax.legend()
    
def isThisQuadCounterclockwise(vertices):
    def crossProduct(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    
    A, B, C, D = vertices
    
    isCounterClockwise = (crossProduct(A, B, C) > 0 and
            crossProduct(B, C, D) > 0 and
            crossProduct(C, D, A) > 0 and
            crossProduct(D, A, B) > 0)
    #print('This quad cell is ordered counterclockwise ', isCounterClockwise)
    return isCounterClockwise

def applyDistanceSizeField(ftLineIds, dx, dxAtBoundary,
                          distMin=5.0, distMax=120.0, refineBoxes=None):
    """Opt-in: express the mesh size law as gmsh fields instead of point sizes.

        h(d) = dx                                          d <= distMin
             = dx + (dxAtBoundary-dx)*(d-distMin)/(distMax-distMin)
             = dxAtBoundary                                d >= distMax

    with d the distance to the nearest fault curve -- linear between two
    clamps (gmsh Threshold), which is what per-point sizing approximates.

    refineBoxes is [(xmin, xmax, ymin, ymax, size), ...], applied as a Min
    against the law, for refining a narrow fault-to-fault gap.

    DEFAULT OFF, and measured as a REGRESSION on xianshuihe.gmsh.lite:

        point sizes (default)   10 triangles,  2 orphans, 34602 cells,  5.7 s
        per-point setSize       14 triangles,  4 orphans, 34866 cells,  ~6 s
        this distance field     24 triangles,  6 orphans, 57249 cells,  3m27s

    Refining a gap does not make Blossom recombination succeed there, it
    multiplies the awkward transitions. Kept because the field form is the
    right way to drive gmsh and may suit other geometries, but do not switch
    a compset to it without re-running checkMeshQuality.

    NOTE any background field overrides per-point sizes, so dxAtBoundary must
    come from the field too, or the whole domain meshes at dx.
    """
    fdist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fdist, "CurvesList", list(ftLineIds))
    gmsh.model.mesh.field.setNumber(fdist, "Sampling", 100)
    fthr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fthr, "InField", fdist)
    gmsh.model.mesh.field.setNumber(fthr, "SizeMin", dx)
    gmsh.model.mesh.field.setNumber(fthr, "SizeMax", dxAtBoundary)
    gmsh.model.mesh.field.setNumber(fthr, "DistMin", distMin)
    gmsh.model.mesh.field.setNumber(fthr, "DistMax", distMax)
    fields = [fthr]
    for (xmin, xmax, ymin, ymax, size) in (refineBoxes or []):
        fb = gmsh.model.mesh.field.add("Box")
        gmsh.model.mesh.field.setNumber(fb, "VIn", size)
        gmsh.model.mesh.field.setNumber(fb, "VOut", dxAtBoundary)
        gmsh.model.mesh.field.setNumber(fb, "XMin", xmin)
        gmsh.model.mesh.field.setNumber(fb, "XMax", xmax)
        gmsh.model.mesh.field.setNumber(fb, "YMin", ymin)
        gmsh.model.mesh.field.setNumber(fb, "YMax", ymax)
        gmsh.model.mesh.field.setNumber(fb, "Thickness", 5.0)
        fields.append(fb)
    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    return fmin


def repairOrphanedSplitNodes(cells, masterSlaveNodeIdRelation, points, ftNodeIds):
    """Give back a cell to any split node left with none.

    replaceMasterWithSlaveNodes swaps master ids to slave ids in every cell
    classified above the fault. Where a node's ENTIRE cell fan is classified
    to one side, the other id ends up in no cell at all -- an orphan, with no
    element to transmit traction.

    R18 removed the commonest cause (side decided by global dy instead of the
    fault normal), but it still happens at the second-to-last node of a fault,
    where the trace ties into a connector line and the local fan is genuinely
    one-sided.

    This re-tests each cell touching the surviving id against the local fault
    tangent and hands the first genuinely-other-side cell back. Returns
    (cells, n_repaired). If no cell qualifies the node is left orphaned rather
    than a cell being moved arbitrarily -- checkMeshQuality will still flag it.
    """
    master, slave = masterSlaveNodeIdRelation
    n_repaired = 0
    for k, (m, s) in enumerate(zip(master, slave)):
        cellNodes = {}
        for iElem, nodeIds in enumerate(cells):
            for nid in nodeIds:
                cellNodes.setdefault(nid, []).append(iElem)
        m_cells = cellNodes.get(m, [])
        s_cells = cellNodes.get(s, [])
        if m_cells and s_cells:
            continue
        if not m_cells and not s_cells:
            continue                      # node is in no cell at all; not ours to fix
        orphan, present, wantAbove = ((m, s, False) if not m_cells else (s, m, True))
        a = max(k - 1, 0)
        b = min(k + 1, len(master) - 1)
        tang = points[ftNodeIds[b], :] - points[ftNodeIds[a], :]
        # Exclude EVERY fault node from the off-fault centroid, not just this
        # one. A cell with an edge on the fault holds two consecutive fault
        # nodes; leaving the second one in drags the centroid along the fault
        # and the side test stops meaning anything.
        ftAll = set(master) | set(slave)
        for iElem in cellNodes.get(present, []):
            onFt = [t for t in cells[iElem] if t in ftAll]
            others = [t for t in cells[iElem] if t not in ftAll]
            if not others or not onFt:
                continue
            quadrant = judgeElemDirect(calcCenterLoc([points[t, :] for t in onFt]),
                                       calcCenterLoc([points[t, :] for t in others]),
                                       ftTangent=tang)
            isAbove = quadrant in (1, 2)
            if isAbove == wantAbove:
                cells[iElem][list(cells[iElem]).index(present)] = orphan
                n_repaired += 1
                break
    return cells, n_repaired


def reorderCellNodesCounterclockwise(cells, pointsWithSplitNodes):
    reorderedCells = cells.copy()

    for cellId, cellNodeIds in enumerate(cells):
        vertices = getCellNodeCoors(cellNodeIds, pointsWithSplitNodes)
        if isThisQuadCounterclockwise(vertices) == False:
            reorderedCells[cellId,0] = cells[cellId,3]
            reorderedCells[cellId,1] = cells[cellId,2]
            reorderedCells[cellId,2] = cells[cellId,1]
            reorderedCells[cellId,3] = cells[cellId,0]
        else:
            noNeedToReorder = True
    return reorderedCells

def calcTanAndLen(xCoors, yCoors):
    # Parametric cubic spline in arc-length s: C^2-smooth tangent, matching
    # Fortran meshgen1.f90's spline-derivative tangent. Returns [tx, ty, len].
    x = np.asarray(xCoors, dtype=float)
    y = np.asarray(yCoors, dtype=float)
    n = len(x)

    seg = np.hypot(np.diff(x), np.diff(y))
    s = np.concatenate(([0.0], np.cumsum(seg)))

    cs_x = CubicSpline(s, x, bc_type="natural")
    cs_y = CubicSpline(s, y, bc_type="natural")
    dxds = cs_x(s, 1)
    dyds = cs_y(s, 1)
    mag = np.hypot(dxds, dyds)
    tx = dxds / mag
    ty = dyds / mag

    # Associated length per node: half of each adjacent segment
    L = np.empty(n)
    L[0]      = 0.5 * seg[0]
    L[-1]     = 0.5 * seg[-1]
    L[1:-1]   = 0.5 * (seg[:-1] + seg[1:])

    return [[float(tx[i]), float(ty[i]), float(L[i])] for i in range(n)]

def writeFilesForEQdyna(pointsWithSplitNodes, cells, masterSlaveNodeIdRelation, ftNodeTanAndLen, ftPhys, modelRange, ftNames):
    
    meshInfo = {}
    meshInfo['totalNumOfNodes'] = len(pointsWithSplitNodes)
    meshInfo['totalNumOfCells'] = len(cells)
    
    for ftName in ftNames:
        meshInfo[ftName] = len(masterSlaveNodeIdRelation[ftName][0])
    
    numOfFtNodes = [] 
    for ftName in meshInfo:
        if ftName in ftNames:
            numOfFtNodes += [meshInfo[ftName]]
    
    maxNumOfFtNodes = max(numOfFtNodes)
    
    nFt = len(ftNames)
    nsmp = np.zeros((maxNumOfFtNodes*nFt,2))
    nsmpTanLen = np.zeros((maxNumOfFtNodes*nFt,3))
    nsmpGeoPhys = np.zeros((maxNumOfFtNodes*nFt,9))
    
    string = ''
    for iFt, ftName in enumerate(ftNames):
        n = len(masterSlaveNodeIdRelation[ftName][0])
        nsmp[iFt*maxNumOfFtNodes:iFt*maxNumOfFtNodes+n, 0:2] = np.array(masterSlaveNodeIdRelation[ftName]).T #+ 1iFt*maxNumOfFtNodes:iFt*maxNumOfFtNodes+n, 0:3    
        nsmpTanLen[iFt*maxNumOfFtNodes:iFt*maxNumOfFtNodes+n, 0:3] = np.array(ftNodeTanAndLen[ftName])
        nsmpGeoPhys[iFt*maxNumOfFtNodes:iFt*maxNumOfFtNodes+n, 0:3] = np.array(ftNodeTanAndLen[ftName])
        nsmpGeoPhys[iFt*maxNumOfFtNodes:iFt*maxNumOfFtNodes+n, 3:9] = np.array(ftPhys[ftName])
        string += str(n)+' ' 
        
    #print(maxNumOfFtNodes)
    #print(meshInfo)
    #print(nsmp)
    
    np.savetxt('vert.txt', pointsWithSplitNodes, fmt='%e')
    np.savetxt('fac.txt', cells, fmt='%d')
    np.savetxt('nsmp.txt', nsmp, fmt='%d')
    np.savetxt('nsmpTanLen.txt', nsmpTanLen, fmt='%e')
    np.savetxt('nsmpGeoPhys.txt', nsmpGeoPhys, fmt='%e')
    
    with open('meshGeneralInfo.txt','w') as f:
        f.write(str(meshInfo['totalNumOfNodes'])+' '+str(meshInfo['totalNumOfCells'])+' \n')
        f.write(string+' \n')
        f.write(str(modelRange['xmin'])+' '+str(modelRange['xmax'])+' '+str(modelRange['ymin'])+' '+str(modelRange['ymax']))
        
    return meshInfo, pointsWithSplitNodes, cells, nsmp, nsmpGeoPhys
