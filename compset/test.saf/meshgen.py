#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!pip install matplotlib meshio
import gmsh
import numpy as np
import meshio
import matplotlib.pyplot as plt
from meshGenLib import *
from userDefinedFaultSysGeoPhys import *

debugMode = False
meshName = 'eqdynaMesh'

# ssaf is split into ssaf1 and ssaf2 for mesh topology.
# In the simulation, they are combined as one fault 'ssaf'.
ftNamesForGmsh = ['ssaf1', 'ssaf2', 'sjfn', 'sjfs']
ftNames = ['ssaf', 'sjfn', 'sjfs']


system = "saf"

# length in km
dx = 0.3
dxAtBoundary = 20
totalSimuTime = 30
vp = 6000
ext = 10 + totalSimuTime*vp/1e3

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
for ftName in ftNamesForGmsh:
    ftFileName = 'user_fault_geometry_input/' + ftName + '.gmt.txt'
    ftLoc = loadFtLoc(ftFileName)
    numOfControlPts, ftEndNodeId, ftRange = createFtPoints(ftLoc, numOfControlPts, dx)
    ftEndNodeIdDict[ftName] = ftEndNodeId
    modelRange = redefineModelRange(modelRange, ftRange)
    lineCount, ftEndLineId = createLinesForFt(ftEndNodeId, lineCount)
    ftEndLineIdDict[ftName] = ftEndLineId
    print(ftLoc)
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

ssaf1Curve = [i+ftEndLineIdDict['ssaf1'][0] for i in range(ftEndLineIdDict['ssaf1'][1]-ftEndLineIdDict['ssaf1'][0]+1)]
ssaf2Curve = [i+ftEndLineIdDict['ssaf2'][0] for i in range(ftEndLineIdDict['ssaf2'][1]-ftEndLineIdDict['ssaf2'][0]+1)]
sjfnCurve = [i+ftEndLineIdDict['sjfn'][0] for i in range(ftEndLineIdDict['sjfn'][1]-ftEndLineIdDict['sjfn'][0]+1)]
sjfsCurve = [i+ftEndLineIdDict['sjfs'][0] for i in range(ftEndLineIdDict['sjfs'][1]-ftEndLineIdDict['sjfs'][0]+1)]
#print(ssaf1Curve, ssaf2Curve, sjfnCurve, sjfsCurve)

# boundary edges
lineCount += 1
T = gmsh.model.geo.addLine(boundaryNodeIdDict['right_top'], boundaryNodeIdDict['left_top'], tag=lineCount)
lineCount += 1
B = gmsh.model.geo.addLine(boundaryNodeIdDict['right_bottom'], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
L = gmsh.model.geo.addLine(boundaryNodeIdDict['left_top'], boundaryNodeIdDict['left_bottom'], tag=lineCount)
lineCount += 1
R = gmsh.model.geo.addLine(boundaryNodeIdDict['right_top'], boundaryNodeIdDict['right_bottom'], tag=lineCount)

# ssaf to boundary
lineCount += 1
S_LT = gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf1'][0], boundaryNodeIdDict['left_top'], tag=lineCount) # ssaf1 start to left top
lineCount += 1
S_LB = gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf1'][0], boundaryNodeIdDict['left_bottom'], tag=lineCount) # ssaf1 start to left bottom
lineCount += 1
S_RT = gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf2'][1], boundaryNodeIdDict['right_top'], tag=lineCount) # ssaf2 end to right top

# ssaf segment connection
lineCount += 1
S1_S2 = gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf1'][1], ftEndNodeIdDict['ssaf2'][0], tag=lineCount) # ssaf1 end to ssaf2 start

# ssaf to sjfn connections
lineCount += 1
S_SJFN = gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf1'][1], ftEndNodeIdDict['sjfn'][0], tag=lineCount) # ssaf1 end to sjfn start
lineCount += 1
S_SJFN_R = gmsh.model.geo.addLine(ftEndNodeIdDict['ssaf2'][1], ftEndNodeIdDict['sjfn'][1], tag=lineCount) # ssaf2 end to sjfn end

# sjfn to sjfs connection
lineCount += 1
SJFN_SJFS = gmsh.model.geo.addLine(ftEndNodeIdDict['sjfn'][1], ftEndNodeIdDict['sjfs'][0], tag=lineCount) # sjfn end to sjfs start

# sjfn/sjfs to boundary
lineCount += 1
SJFS_RB = gmsh.model.geo.addLine(ftEndNodeIdDict['sjfs'][1], boundaryNodeIdDict['right_bottom'], tag=lineCount) # sjfs end to right bottom
lineCount += 1
SJFN_LB = gmsh.model.geo.addLine(ftEndNodeIdDict['sjfn'][0], boundaryNodeIdDict['left_bottom'], tag=lineCount) # sjfn start to left bottom

# ssaf-top block
surfaceCount = createSurface(ssaf1Curve+[S1_S2]+ssaf2Curve+[S_RT, T, -S_LT], surfaceCount)
# left block
surfaceCount = createSurface([S_LB, -L, -S_LT], surfaceCount)
# ssaf-sjfn wedge block
surfaceCount = createSurface([S1_S2]+ssaf2Curve+[S_SJFN_R]+addMinusToList(sjfnCurve[::-1])+[-S_SJFN], surfaceCount)
# bottom-left block
surfaceCount = createSurface(addMinusToList(ssaf1Curve[::-1])+[S_LB, -SJFN_LB, -S_SJFN], surfaceCount)
# sjf-bottom block
surfaceCount = createSurface([SJFN_LB, -B, -SJFS_RB]+addMinusToList(sjfsCurve[::-1])+[-SJFN_SJFS]+addMinusToList(sjfnCurve[::-1]), surfaceCount)
# right block
surfaceCount = createSurface([S_SJFN_R, SJFN_SJFS]+sjfsCurve+[SJFS_RB, -R, -S_RT], surfaceCount)

gmsh.model.geo.synchronize()

for iSur in range(surfaceCount):
    surfaceId = iSur+1
    print('recombining surface id ', surfaceId)
    gmsh.model.mesh.setRecombine(2, surfaceId)

#gmsh.model.mesh.coherence() not available in python
gmsh.model.mesh.generate(2)
#gmsh.model.mesh.removeUnusedEntities()
gmsh.write('fem_mesh_output/' + meshName+'.msh')
gmsh.option.setNumber("Mesh.Smoothing", 2)
#gmsh.finalize()


# In[3]:


ftTag = {}
ftTag['ssaf'] = ssaf1Curve+[S1_S2]+ssaf2Curve
ftTag['sjfn'] = sjfnCurve
ftTag['sjfs'] = sjfsCurve
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

fig, ax = plt.subplots(figsize=(15, 20), dpi=600)
ax.scatter(points[:, 0], points[:, 1], s=0.01, color='red', zorder=1)
for cell in cells:
    vertices = points[cell]
    ax.add_patch(plt.Polygon(vertices, edgecolor='black', linewidth=0.1, fill=False))
ax.set_aspect('equal')

plt.savefig('fem_mesh_output/meshWOSplitNode.png', dpi=600)
if debugMode==True:
    plt.show()


# In[6]:


ftNodeIdsDict={}
for key in ftNames:
    ftNodeIdsDict[key] = locateFtNodeIds(points, xCoorDict[key], yCoorDict[key], tolerance)
    if debugMode==True:
        print(ftNodeIdsDict[key])


# In[7]:


# check and plot ft nodes
fig, ax = plt.subplots(figsize=(15, 10), dpi=600)

ax.scatter(points[:, 0], points[:, 1], s=0.1, color='red', zorder=1)
for ftName in ftNames:
    ax.scatter(points[ftNodeIdsDict[ftName],0], points[ftNodeIdsDict[ftName],1], s=0.3, color='black', zorder=2)

plt.savefig('fem_mesh_output/meshWithFaultNodes.png', dpi=600)
if debugMode==True:
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


# check and plot ft nodes
fig, ax = plt.subplots(figsize=(15, 10), dpi=600)
plt.scatter(points[:, 0], points[:, 1], s=0.01, color='red', zorder=1)
for key in ftNames:
    plt.scatter(pointsWithSplitNodes[ftNodeIdsDict[key],0], pointsWithSplitNodes[ftNodeIdsDict[key],1], s=0.03, color='black', zorder=2)
    plt.scatter(pointsWithSplitNodes[slaveNodeIdsDict[key],0], pointsWithSplitNodes[slaveNodeIdsDict[key],1], marker='*', s=0.003, color='blue', zorder=2)
#plt.savefig('mesh.png', dpi=600)
if debugMode==True:
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


ftNodeTanAndLen = {}
for ftNameKey in ftNames:
    ftNodeTanAndLen[ftNameKey] = calcTanAndLen(xCoorDict[ftNameKey], yCoorDict[ftNameKey])


# In[16]:


ftPhys = defineSysPhys(system, ftNames, xCoorDict, yCoorDict)


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


if debugMode==True:
    plotSystemPhys(ftPhys, ftNames, xCoorDict)
