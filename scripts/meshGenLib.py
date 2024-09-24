#! /usr/env/bin python3
import numpy as np
import gmsh
import matplotlib.pyplot as plt

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

def judgeElemDirect(cxOnFt, cxOffFt):
    quadrant = 0
    elemDirectVec = np.array(cxOffFt) - np.array(cxOnFt)
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
    elemIdsAboveFt = []
    elemIdsBelowFt = []
    for iNode in range(len(ftNodeIds)-1):
        pairOfNodes = [ftNodeIds[iNode], ftNodeIds[iNode+1]]
        for iElem, elemTag in enumerate(cells):
            nodeIds = cells[iElem]
            if pairOfNodes[0] in nodeIds and pairOfNodes[1] in nodeIds:
                otherTwoNodes = [tag for tag in nodeIds if tag not in pairOfNodes]
                #print(nodeIds, pairOfNodes, otherTwoNodes)
                
                coorsOnFt = [points[nodeId,:] for nodeId in pairOfNodes]
                coorsOffFt = [points[nodeId,:] for nodeId in otherTwoNodes]
                cxOnFt = calcCenterLoc(coorsOnFt)
                cxOffFt = calcCenterLoc(coorsOffFt)
                #print(cxOnFt, cxOffFt)
                quadrant = judgeElemDirect(cxOnFt, cxOffFt)
                #print(quadrant)

                if quadrant==1 or quadrant==2:
                    elemIdsAboveFt += [iElem]  
                    #gmsh.model.addPhysicalGroup(2, [elemTag], physicalGroupForElemAboveFt)
                    #gmsh.model.setPhysicalName(2, physicalGroupForElemAboveFt, 'elemAboveFt')

                if quadrant==3 or quadrant==4:
                    elemIdsBelowFt += [iElem]
                    #gmsh.model.addPhysicalGroup(2, [elemTag], physicalGroupForElemBelowFt)
                    #gmsh.model.setPhysicalName(2, physicalGroupForElemBelowFt, 'elemBelowFt')
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
    def calcLen(x1,y1, x2,y2):
        return ((x2-x1)**2+(y2-y1)**2)**0.5
    
    tan = []
    # for i = 0 
    ftNodeLength =  0.5*calcLen(xCoors[1],yCoors[1],xCoors[0],yCoors[0])
    tangent = [(xCoors[1]-xCoors[0])/ftNodeLength/2, (yCoors[1]-yCoors[0])/ftNodeLength/2, ftNodeLength]
    tan += [tangent]
    
    for i in range(1, len(xCoors) - 1):
        # Calculate the tangent vector as the difference between consecutive points
        ftNodeLength = 0.5*(calcLen(xCoors[i+1],yCoors[i+1],xCoors[i],yCoors[i]) + calcLen(xCoors[i-1],yCoors[i-1],xCoors[i],yCoors[i]))
        
        len1 = calcLen(xCoors[i+1], yCoors[i+1], xCoors[i-1], yCoors[i-1])
        tangent = [(xCoors[i+1]-xCoors[i-1])/len1, (yCoors[i+1]-yCoors[i-1])/len1, ftNodeLength]
        tan += [tangent]
        
    # for i = len(xCoors)
    ftNodeLength =  0.5*calcLen(xCoors[len(xCoors)-1],yCoors[len(xCoors)-1],xCoors[len(xCoors)-2],yCoors[len(xCoors)-2])
    tangent = [(xCoors[len(xCoors)-1]-xCoors[len(xCoors)-2])/ftNodeLength/2, (yCoors[len(xCoors)-1]-yCoors[len(xCoors)-2])/ftNodeLength/2, ftNodeLength]
    tan += [tangent]
        
        ## Rotate the tangent vector by 90 degrees to get the normal vector
        #normal = np.array([-tangent[1], tangent[0]])
        # Normalize the normal vector
        #normal = normal / np.linalg.norm(normal)
        #normals.append(normal)
    return tan

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
    
    nsmp = np.zeros((maxNumOfFtNodes*3,2))
    nsmpTanLen = np.zeros((maxNumOfFtNodes*3,3))
    nsmpGeoPhys = np.zeros((maxNumOfFtNodes*3,9))
    
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
    np.savetxt('nsmpTanlen.txt', nsmpTanLen, fmt='%e')
    np.savetxt('nsmpGeoPhys.txt', nsmpGeoPhys, fmt='%e')
    
    with open('meshGeneralInfo.txt','w') as f:
        f.write(str(meshInfo['totalNumOfNodes'])+' '+str(meshInfo['totalNumOfCells'])+' \n')
        f.write(string+' \n')
        f.write(str(modelRange['xmin'])+' '+str(modelRange['xmax'])+' '+str(modelRange['ymin'])+' '+str(modelRange['ymax']))
        
    return meshInfo, pointsWithSplitNodes, cells, nsmp, nsmpGeoPhys

def defineSysPhys(ftSystem, ftNames, xCoorDict, yCoorDict):
    ftPhys = {}
    for ftNameKey in ftNames:
        tmp = []
        for i, xcoor in enumerate(xCoorDict[ftNameKey]):
            
            if ftSystem=="none":
                ftType = 1 # 1: left-strike; -1: right-strike; 2: thrust; -2: normal. 
                ftDip = 90 # 90: strike-slip; positive: tilting to y+; negative; tilting to y-.
                ftLoadMaxShear = 1.43e-14
                ftLoadAngle = -999
                ftLoadWt = 1.
                ftVis = 6e21 
                tmp += [[ftType, ftDip, ftLoadMaxShear, ftLoadAngle, ftLoadWt, ftVis]]
            elif ftSystem=="subei":
                ftType = 1
                ftDip = 90
                ftLoadMaxShear = 1.43e-14
                ftLoadAngle = -999
                ftLoadWt = 1.
                ftVis = 6e21
                
                if ftNameKey=='sbt':
                    ftType = 2   
                    ftDip = 30
                
                if ftNameKey=='atf' and xcoor>0:
                    ftLoadWt = .8
                
                if ftNameKey=='dxs':
                    ftLoadWt = .1
                
                tmp += [[ftType, ftDip, ftLoadMaxShear, ftLoadAngle, ftLoadWt, ftVis]]
                
        ftPhys[ftNameKey] = tmp
    return ftPhys

def plotSystemPhys(ftPhys, ftNames, xCoorDict):
    fig, ax = plt.subplots(3, 2, figsize=(15, 10), dpi=600)
    
    for key in ftNames:
        ax[0,0].plot(xCoorDict[key], [row[0] for row in ftPhys[key]]) 
        ax[0,1].plot(xCoorDict[key], [row[1] for row in ftPhys[key]])         
        ax[1,0].plot(xCoorDict[key], [row[2] for row in ftPhys[key]]) 
        ax[1,1].plot(xCoorDict[key], [row[3] for row in ftPhys[key]]) 
        ax[2,0].plot(xCoorDict[key], [row[4] for row in ftPhys[key]]) 
        ax[2,1].plot(xCoorDict[key], [row[5] for row in ftPhys[key]])         
    #plt.savefig('mesh.png', dpi=600)
    plt.show()