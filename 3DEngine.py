#Thien Vu
#import stuff
import math
import copy
import sys
from tkinter import *
#canvas sizes
CanvasWidth = 600
CanvasHeight = 600
#distance of viewpoint
d = 500
#ambient and pointlight intensity
Ia = 0.3
Ip = 0.7
#object intensity for every object
Kd = 0.6
Ks = 1-Kd
#lighting vector and view vector
L = [1,1,-1]
V = [0,0,-1]

#specular index
specIndex = 8

#edge class
#used in edge table
class Edge():
    xStart = 0
    yStart = 0
    yEnd = 0
    dX = 0
    zStart = 0
    dZ = 0
    xEnd = 0
    iStart = []
    iEnd = []
    dI = []
    nStart = []
    nEnd = []
    dN = []
    #construction function
    def __init__(self, xStart, yStart, yEnd, dX, zStart, dZ, xEnd, zEnd):
        self.xStart = xStart
        self.yStart = yStart
        self.yEnd = yEnd
        self.dX = dX
        self.zStart = zStart
        self.dZ = dZ
        self.xEnd = xEnd
        self.zEnd = zEnd

    #set intensities for gouraud shading
    def setIntensities(self, iStart, iEnd, dI):
        self.iStart = iStart
        self.iEnd = iEnd
        self.dI = dI

    # set normals for phong shading
    def setNormals(self, nStart, nEnd, dN):
        self.nStart = nStart
        self.nEnd = nEnd
        self.dN = dN
# ***************************** Initialize Pyramid Object ***************************
# Definition  of the five underlying points
apex = [(-200+200),50,100-150]
base1 = [(-150+200),-50,50-150]
base2 = [(-150+200),-50,150-150]
base3 = [(-250+200),-50,150-150]
base4 = [(-250+200),-50,50-150]


# Definition of the five polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
frontpoly = [apex,base1,base4]
rightpoly = [apex,base2,base1]
backpoly = [apex,base3,base2]
leftpoly = [apex,base4,base3]
bottompoly = [base1,base2,base3,base4]

# Definition of the object
Pyramid = [bottompoly, frontpoly, rightpoly, backpoly, leftpoly]
# Definition of the Pyramid's underlying point cloud.  No structure, just the points.
PyramidPointCloud = [apex, base1, base2, base3, base4]
DefaultPyramidPointCloud = copy.deepcopy(PyramidPointCloud)

# Definition of the color associated with each face of the pyramid
PyramidColor = ["black", "red", "green", "blue", "yellow"]
#************************************************************************************

# ***************************** Initialize Cube Object ***************************
# Definition  of the 8 underlying points
frontTRCube = [(250-200),100-50,0]
frontTLCube = [(150-200),100-50,0]
frontBRCube = [(250-200),0-50,0]
frontBLCube = [(150-200),0-50,0]
backTRCube = [(250-200),100-50,100]
backTLCube = [(150-200),100-50,100]
backBRCube = [(250-200),0-50,100]
backBLCube = [(150-200),0-50,100]


# Definition of the 6 polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
frontCube = [frontTLCube, frontTRCube, frontBRCube, frontBLCube]
backCube = [backTRCube, backTLCube, backBLCube, backBRCube]
leftCube = [backTLCube, frontTLCube, frontBLCube, backBLCube]
rightCube = [frontTRCube, backTRCube, backBRCube, frontBRCube]
topCube = [backTLCube, backTRCube, frontTRCube, frontTLCube]
bottomCube = [frontBLCube, frontBRCube, backBRCube, backBLCube]

# Definition of the object
Cube = [frontCube, backCube, leftCube, rightCube, topCube, bottomCube]

# Definition of the Cube's underlying point cloud.  No structure, just the points.
CubePointCloud = [frontTRCube, frontTLCube, frontBRCube, frontBLCube, backTRCube, backTLCube, backBRCube, backBLCube]
DefaultCubePointCloud = copy.deepcopy(CubePointCloud)
CubeColor = ["white", "#cccccc", "#999999", "#666666", "#333333", "black"]
#************************************************************************************

# ***************************** Initialize Cube2 Object ***************************
# Definition  of the 8 underlying points
frontTRCube2 = [75,75,-75]
frontTLCube2 = [-75,75,-75]
frontBRCube2 = [75,-75,-75]
frontBLCube2 = [-75,-75,-75]
backTRCube2 = [75,75,75]
backTLCube2 = [-75,75,75]
backBRCube2 = [75,-75,75]
backBLCube2 = [-75,-75,75]


# Definition of the 6 polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
frontCube2 = [frontTLCube2, frontTRCube2, frontBRCube2, frontBLCube2]
backCube2 = [backTRCube2, backTLCube2, backBLCube2, backBRCube2]
leftCube2 = [backTLCube2, frontTLCube2, frontBLCube2, backBLCube2]
rightCube2 = [frontTRCube2, backTRCube2, backBRCube2, frontBRCube2]
topCube2 = [backTLCube2, backTRCube2, frontTRCube2, frontTLCube2]
bottomCube2 = [frontBLCube2, frontBRCube2, backBRCube2, backBLCube2]

# Definition of the object
Cube2 = [frontCube2, backCube2, leftCube2, rightCube2, topCube2, bottomCube2]

# Definition of the Cube's underlying point cloud.  No structure, just the points.
Cube2PointCloud = [frontTRCube2, frontTLCube2, frontBRCube2, frontBLCube2, backTRCube2, backTLCube2, backBRCube2, backBLCube2]
DefaultCube2PointCloud = copy.deepcopy(Cube2PointCloud)
#color
Cube2Color = ["white", "#cccccc", "#999999", "#666666", "#333333", "black"]
#************************************************************************************

# ***************************** Initialize HexagonalPrism Object ***************************
# Definition  of the 12 underlying points
frontTLHexP = [-20, 135, 0]
frontTRHexP = [20, 135, 0]
frontMRHexP = [40, 100, 0]
frontBRHexP = [20, 65, 0]
frontBLHexP = [-20, 65, 0]
frontMLHexP = [-40, 100, 0]
backTLHexP = [-20, 135, 50]
backTRHexP = [20, 135, 50]
backMRHexP = [40, 100, 50]
backBRHexP = [20, 65, 50]
backBLHexP = [-20, 65, 50]
backMLHexP = [-40, 100, 50]


# Definition of the 8 polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
frontHexP = [frontTLHexP, frontTRHexP, frontMRHexP, frontBRHexP, frontBLHexP, frontMLHexP]
backHexP = [backMLHexP, backBLHexP, backBRHexP, backMRHexP, backTRHexP, backTLHexP]
topLHexP = [backTLHexP, frontTLHexP, frontMLHexP, backMLHexP]
topHexP = [frontTLHexP, backTLHexP, backTRHexP, frontTRHexP]
topRHexP = [frontTRHexP, backTRHexP, backMRHexP, frontMRHexP]
botRHexP = [frontMRHexP, backMRHexP, backBRHexP, frontBRHexP]
botHexP = [frontBRHexP, backBRHexP, backBLHexP, frontBLHexP]
botLHexP = [backMLHexP, frontMLHexP, frontBLHexP, backBLHexP]


# Definition of the object
HexP = [frontHexP, backHexP, topLHexP, topHexP, topRHexP, botRHexP, botHexP, botLHexP]

# Definition of the Cube's underlying point cloud.  No structure, just the points.
HexPPointCloud = [frontTLHexP, frontTRHexP, frontMRHexP, frontBRHexP, frontBLHexP, frontMLHexP,backTLHexP, backTRHexP, backMRHexP, backBRHexP, backBLHexP, backMLHexP]
DefaultHexPPointCloud = copy.deepcopy(HexPPointCloud)
HexPColor = ["RoyalBlue1","RoyalBlue2","DodgerBlue2","DodgerBlue3","DodgerBlue4","SteelBlue1","SteelBlue2","SteelBlue3"]
#************************************************************************************

# ***************************** Initialize Octo-Cylinder Object ***************************
# Definition of the 16 underlying points
octoCylFront1 = [-50,120.7107,50]
octoCylFront2 = [50,120.7107,50]
octoCylFront3 = [120.7107,50,50]
octoCylFront4 = [120.7107,-50,50]
octoCylFront5 = [50,-120.7107,50]
octoCylFront6 = [-50,-120.7107,50]
octoCylFront7 = [-120.7107,-50,50]
octoCylFront8 = [-120.7107,50,50]
octoCylBack1 = [-50,120.7107,450]
octoCylBack2 = [50,120.7107,450]
octoCylBack3 = [120.7107,50,450]
octoCylBack4 = [120.7107,-50,450]
octoCylBack5 = [50,-120.7107,450]
octoCylBack6 = [-50,-120.7107,450]
octoCylBack7 = [-120.7107,-50,450]
octoCylBack8 = [-120.7107,50,450]
# Definition of the ten polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
octoCylNorthPoly = [octoCylFront1, octoCylBack1, octoCylBack2, octoCylFront2]
octoCylNorthEastPoly = [octoCylFront2, octoCylBack2, octoCylBack3, octoCylFront3]
octoCylEastPoly = [octoCylFront3, octoCylBack3, octoCylBack4, octoCylFront4]
octoCylSouthEastPoly = [octoCylFront4, octoCylBack4, octoCylBack5, octoCylFront5]
octoCylSouthPoly = [octoCylFront5, octoCylBack5, octoCylBack6, octoCylFront6]
octoCylSouthWestPoly = [octoCylFront6, octoCylBack6, octoCylBack7, octoCylFront7]
octoCylWestPoly = [octoCylFront7, octoCylBack7, octoCylBack8, octoCylFront8]
octoCylNorthWestPoly = [octoCylFront8, octoCylBack8, octoCylBack1, octoCylFront1]
octoCylFrontPoly = [octoCylFront1, octoCylFront2, octoCylFront3, octoCylFront4, octoCylFront5, octoCylFront6, octoCylFront7, octoCylFront8]
octoCylBackPoly = [octoCylBack1, octoCylBack8, octoCylBack7, octoCylBack6, octoCylBack5, octoCylBack4, octoCylBack3, octoCylBack2]
# Definition of the cylinder object
cylinder = [octoCylNorthPoly, octoCylNorthEastPoly, octoCylEastPoly, octoCylSouthEastPoly, octoCylSouthPoly, octoCylSouthWestPoly, octoCylWestPoly, octoCylNorthWestPoly,
octoCylFrontPoly, octoCylBackPoly]
OctoCylPointCloud = [octoCylFront1, octoCylFront2, octoCylFront3, octoCylFront4, octoCylFront5, octoCylFront6, octoCylFront7, octoCylFront8,
                     octoCylBack1, octoCylBack2, octoCylBack3, octoCylBack4, octoCylBack5, octoCylBack6, octoCylBack7, octoCylBack8]
DefaultOctoCylPointCloud = copy.deepcopy(OctoCylPointCloud)
OctoCylColor = ["#006400","#006400","#006400","#006400","#006400","#006400","#006400","#006400","#006400","#006400"]
#list of objects
#objectList = [Pyramid, Cube, Cube2, HexP]
objectList = [cylinder]
#objectPointCloudList = [PyramidPointCloud, CubePointCloud, Cube2PointCloud, HexPPointCloud]
objectPointCloudList = [OctoCylPointCloud]
#objectDefaultPointCloudList = [DefaultPyramidPointCloud, DefaultCubePointCloud, DefaultCube2PointCloud,DefaultHexPPointCloud]
objectDefaultPointCloudList = [DefaultOctoCylPointCloud]
#objectColorList = [PyramidColor, CubeColor, Cube2Color, HexPColor]
objectColorList = [OctoCylColor]
currentObj = [0]

#zbuffer setup
zBuffer = [[d for i in range(CanvasHeight)] for j in range(CanvasWidth)]
defaultZBuffer = copy.deepcopy(zBuffer)

#render mode
renderM = [1]

# This function resets the pyramid to its original size and location in 3D space
# Note that you have to be careful to update the values in the existing PyramidPointCloud
# structure rather than creating a new structure or just switching a pointer.  In other
# words, you'll need manually update the value of every x, y, and z of every point in
# point cloud (vertex list).
def resetPyramid():
    for i in range(len(objectDefaultPointCloudList[currentObj[0]])):
        for j in range(3):
            objectPointCloudList[currentObj[0]][i][j] = objectDefaultPointCloudList[currentObj[0]][i][j]


# This function translates an object by some displacement.  The displacement is a 3D
# vector so the amount of displacement in each dimension can vary.
def translate(object, displacement):
    # for each point
    for i in range(len(object)):
        # for x,y,z
        for j in range(len(object[i])):
            object[i][j] += displacement[j]

# This function performs a simple uniform scale of an object assuming the object is
# centered at the origin.  The scalefactor is a scalar.
def scale(object, scalefactor):
    refP, negRefP = getRef(object)
    translate(object, negRefP)
    # for each point
    for i in range(len(object)):
        # for x,y,z
        for j in range(len(object[i])):
            object[i][j] *= scalefactor
    translate(object, refP)


# This function performs a rotation of an object about the Z axis (from +X to +Y)
# by 'degrees', assuming the object is centered at the origin.  The rotation is CCW
# in a LHS when viewed from -Z [the location of the viewer in the standard position]
def rotateZ(object, degrees):
    refP, negRefP = getRef(object)
    translate(object, negRefP)

    # convert degrees to radians
    rad = degrees * math.pi / 180
    # for each point
    for i in range(len(object)):
        # keep original coords
        coords = []
        for j in range(len(object[i])):
            coords.append(object[i][j])
        # change x and y, z stays the same
        object[i][0] = coords[0] * math.cos(rad) - coords[1] * math.sin(rad)
        object[i][1] = coords[0] * math.sin(rad) + coords[1] * math.cos(rad)
        # dont change z
    translate(object, refP)


# This function performs a rotation of an object about the Y axis (from +Z to +X)
# by 'degrees', assuming the object is centered at the origin.  The rotation is CW
# in a LHS when viewed from +Y looking toward the origin.
def rotateY(object, degrees):
    refP, negRefP = getRef(object)
    translate(object, negRefP)
    # conv deg to rad
    rad = degrees * math.pi / 180
    for i in range(len(object)):
        # keep original coords
        coords = []
        for j in range(len(object[i])):
            coords.append(object[i][j])
        # change x and z, y stays the same
        object[i][0] = coords[0] * math.cos(rad) + coords[2] * math.sin(rad)
        # dont change y
        object[i][2] = -coords[0] * math.sin(rad) + coords[2] * math.cos(rad)
    translate(object, refP)


# This function performs a rotation of an object about the X axis (from +Y to +Z)
# by 'degrees', assuming the object is centered at the origin.  The rotation is CW
# in a LHS when viewed from +X looking toward the origin.
def rotateX(object, degrees):
    refP, negRefP = getRef(object)
    translate(object, negRefP)
    # conv deg to rad
    rad = degrees * math.pi / 180
    for i in range(len(object)):
        # keep original coords
        coords = []
        for j in range(len(object[i])):
            coords.append(object[i][j])
        # change y and z, x stays the same
        # dont change x
        object[i][1] = coords[1] * math.cos(rad) - coords[2] * math.sin(rad)
        object[i][2] = coords[1] * math.sin(rad) + coords[2] * math.cos(rad)
    translate(object, refP)

#1 is wireframe
#2 is polyfill
# The function will draw an object by repeatedly calling drawPoly on each polygon in the object
def drawObject(object, color, pColor, zBuffer):
    #draw each polygon for the object
    objNorms = []
    newONorm = []
    visList = []
    #get norms and visibility of each polygon
    for i in range(len(object)):
        vis, norm = cullBF([0,0,-d], object[i])
        objNorms.append(norm)
        visList.append(vis)
    if renderM[0] > 4:
        #this is a cylinder-specific way of computing new norms
        #for each poly(except caps)
        for i in range(8):
            newPNorm = []
            #get index of previous and next poly
            pPoly = i-1
            if pPoly < 0:
                pPoly = 7
            nPoly = (i+1)%8
            #for each vertex
            for j in range(4):
                #for each dimension
                VNorm = []
                for k in range(3):
                    #combine the vectors
                    if j < 2:
                        VNorm.append(objNorms[i][j][k] + objNorms[pPoly][3-j][k])
                    else:
                        VNorm.append(objNorms[i][j][k] + objNorms[nPoly][3-j][k])
                #normalize the vector
                newPNorm.append(normalize(VNorm))
            newONorm.append(newPNorm)
            #print(f"Polygon {i}: {newONorm[i]}")
        #add the last two normals of the caps
        newONorm.append(objNorms[8])
        newONorm.append(objNorms[9])

    #draw poly if not culled
    for i in range(len(visList)):
        if visList[i]:
            if renderM[0] < 5:
                drawPoly(object[i], color, pColor[i], zBuffer, objNorms[i])
            else:
                drawPoly(object[i], color, pColor[i], zBuffer, newONorm[i])

# Draws every object in objectList
# selected object will be drawn in red
def drawObjects(objectList, objectColorList, zBuffer):
    resetZBuffer()
    #for each object in object list
    for i in range(len(objectList)):
        if i == currentObj[0]:
            drawObject(objectList[i], 'black', objectColorList[i], zBuffer)
        else:
            drawObject(objectList[i], 'black', objectColorList[i], zBuffer)
#resets zbuffer
#used to when drawing objects
def resetZBuffer():
    for i in range(len(zBuffer)):
        for j in range(len(zBuffer[0])):
            zBuffer[i][j] = defaultZBuffer[i][j]

#selects the next or previous object in list
#self explanatory tbh
def changeObj(num):
    #modify object position
    currentObj[0] += num
    #if it goes over the amount of objects, loop back to zero
    if currentObj[0] > len(objectList) - 1:
        currentObj[0] = 0
    #if it goes back lower than zero, go to last object in list
    if currentObj[0] < 0:
        currentObj[0] = len(objectList) - 1

# This function will draw a polygon by repeatedly calling drawLine on each pair of points
# making up the object.  Remember to draw a line between the last point and the first.
def drawPoly(poly, color, pColor, zBuffer, polyNorm):
    #initiate display edge array
    dEdges = []
    #fill edge array
    #if wireframe
    if renderM[0] == 1:
        # draw each line for the polygon
        for i in range(len(poly) - 1):
            drawLine(poly[i], poly[i + 1], color)
        drawLine(poly[len(poly) - 1], poly[0], color)
    elif renderM[0] > 4:
        for i in range(len(poly) - 1):
            dEdges.append(getEdge(poly[i], poly[i+1], color, polyNorm[i], polyNorm[i+1]))
        dEdges.append(getEdge(poly[len(poly) - 1], poly[0], color, polyNorm[len(poly) - 1], polyNorm[0]))
    else:
        for i in range(len(poly) - 1):
            dEdges.append(getEdge(poly[i], poly[i+1], color, 0, 0))
        dEdges.append(getEdge(poly[len(poly) - 1], poly[0], color, 0, 0))
    polyFill(dEdges, pColor, zBuffer, color, polyNorm)

#calculates the normal(s) for a polygon
def calcNorms(verts):
    norms = []
    vLen = len(verts)
    #for each vertex in a polygon
    for i in range(vLen):
        # get p and q, which are vectors from vertice 0 to vertice 1, and vertice 0 to vertice 2
        p = [verts[(i+1)%vLen][0] - verts[i][0], verts[(i+1)%vLen][1] - verts[i][1], verts[(i+1)%vLen][2] - verts[i][2]]
        q = [verts[(i+2)%vLen][0] - verts[i][0], verts[(i+2)%vLen][1] - verts[i][1], verts[(i+2)%vLen][2] - verts[i][2]]
        # get surface normal of polygon
        n = [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]]

        # normalize the vector to a unit vector
        uV = normalize(n)
        norms.append(uV)
        #if its not gouraud or phong shading, only 1 normal is needed for backface culling
        if renderM[0] < 5:
            break
    return norms

#culls backfacing polygons from viewpoint
#takes in the viewpoint, which is a coordinate
#and verts, which is a list of 3 coordinates
def cullBF(viewpoint, verts):
    #calculate the norms
    norms = calcNorms(verts)
    #get first normal for backface culling
    uV = norms[0]

    # get plane offset
    pOffset = uV[0] * verts[0][0] + uV[1] * verts[0][1] + uV[2] * verts[0][2]

    #finally, check if visible
    if uV[0] * viewpoint[0] + uV[1] * viewpoint[1] + uV[2] * viewpoint[2] - pOffset > 0:
        return True, norms
    else:
        return False, norms

#gets the edge table
#also switches start and end coordinates if needed
#and sorts edges before making the edge table
def calcEdgeTable(dEdges):
    #switch start and end coordinates if start Y > end Y
    for edge in dEdges:
        # if the y in start edge > y in end edge
        if edge[0][1] > edge[1][1]:
            # switch the whole coordinate
            edge[0], edge[1] = edge[1], edge[0]
    # sort edges by increasing y start values
    for i in range(len(dEdges)):
        for j in range(0, len(dEdges) - i - 1):
            if dEdges[j][0][1] > dEdges[j + 1][0][1]:
                dEdges[j], dEdges[j + 1] = dEdges[j + 1], dEdges[j]
    # make the edge table
    edgeT = []
    # fill with edge table info
    for edge in dEdges:
        if edge[0][1] != edge[1][1]:
            # get all the properties of the edge table
            xStart = edge[0][0]
            yStart = edge[0][1]
            yEnd = edge[1][1]
            dX = (edge[1][0] - edge[0][0]) / (edge[1][1] - edge[0][1])
            zStart = edge[0][2]
            dZ =(edge[1][2] - edge[0][2]) / (edge[1][1] - edge[0][1])
            xEnd = edge[1][0]
            zEnd = edge[1][2]
            edgeInfo = Edge(xStart, yStart, yEnd, dX, zStart, dZ, xEnd, zEnd)
            #calculate dIntensities for gouraud shading
            if renderM[0] == 5:
                intensityStart = calcAmbDiffSpec(Ia, Kd, Ks, Ip, edge[0][3], L, 1, V, specIndex)
                intensityEnd = calcAmbDiffSpec(Ia, Kd, Ks, Ip, edge[1][3], L, 1, V, specIndex)
                dIX = (intensityEnd[0] - intensityStart[0])/ (yEnd - yStart)
                dIY = (intensityEnd[1] - intensityStart[1])/ (yEnd - yStart)
                dIZ = (intensityEnd[2] - intensityStart[2])/ (yEnd - yStart)
                edgeInfo.setIntensities(intensityStart, intensityEnd, [dIX, dIY, dIZ])
                #print(f"Starting intensities: {intensityStart} Ending intensities: {intensityEnd}")
            #calculate dNormals for phong shading
            if renderM[0] == 6:
                nStart = edge[0][3]
                nEnd = edge[1][3]
                dNX = (nEnd[0] - nStart[0]) / (yEnd - yStart)
                dNY = (nEnd[1] - nStart[1]) / (yEnd - yStart)
                dNZ = (nEnd[2] - nStart[2]) / (yEnd - yStart)
                edgeInfo.setNormals(normalize(nStart), normalize(nEnd), [dNX, dNY, dNZ])
            edgeT.append(edgeInfo)
    return edgeT

#fills the coordinate polygon
#takes edge list and polygon color as inputs
def polyFill(dEdges, pColor, zBuffer, color, norm):
    #comput flat shading color
    if renderM[0] == 4:
        flatColor = triColorHexCode(calcAmbDiffSpec(Ia, Kd, Ks, Ip, normalize(norm[0]), L, 1, V, specIndex))
    edgeThickness = 2
    #make edge table with xStart, yStart, yEnd, dX, zStart, and dZ
    edgeT = calcEdgeTable(dEdges)
    #if table is empty polygon is less than a pixel long
    if len(edgeT) == 0:
        return

    #get vertical boundaries
    topY = edgeT[0].yStart
    botY = edgeT[-1].yEnd
    #get horizontal boundaries

    #initialize indices of edges being worked on and next edge in line
    i = 0
    j = 1
    k = 2

    #initial x,z, and intensity boundaries
    eIX = edgeT[i].xStart
    eJX = edgeT[j].xStart
    eIZ = edgeT[i].zStart
    eJZ = edgeT[j].zStart
    eIXEnd = edgeT[i].xEnd
    eJXEnd = edgeT[j].xEnd
    #boundaries for gouraud and phong shading
    if renderM[0] == 5:
        eII = []
        eJI = []
        for c in range(3):
            eII.append(edgeT[i].iStart[c])
            eJI.append(edgeT[j].iStart[c])
    if renderM[0] == 6:
        eIN = []
        eJN = []
        for c in range(3):
            eIN.append(edgeT[i].nStart[c])
            eJN.append(edgeT[j].nStart[c])

    #while polygon isnt filled all the way at the bottom
    y = topY
    while y < botY + 1:
        #get left and right edges
        if eIX < eJX:
            leftX = eIX
            rightX = eJX
            leftZ = eIZ
            rightZ = eJZ
            leftBound = eIXEnd
            rightBound = eJXEnd
            if renderM[0] == 5:
                leftI = eII
                rightI = eJI
            if renderM[0] == 6:
                leftN = eIN
                rightN = eJN
        else:
            leftX = eJX
            rightX = eIX
            leftZ = eJZ
            rightZ = eIZ
            leftBound = eJXEnd
            rightBound = eIXEnd
            if renderM[0] == 5:
                leftI = eJI
                rightI = eII
            if renderM[0] == 6:
                leftN = eJN
                rightN = eIN
        #calculate dZ for this y level
        #also dI and dN for gouraud and phong if needed
        z = leftZ
        if renderM[0] == 5:
            intsy = []
            for c in range(len(leftI)):
                intsy.append(leftI[c])
        if renderM[0] == 6:
            normsy = []
            for c in range(len(leftN)):
                normsy.append(leftN[c])
        if rightX - leftX != 0:
            dZCurr = (float(rightZ) - leftZ) / (rightX - leftX)
            if renderM[0] == 5:
                dICurr = []
                for c in range(3):
                    dICurr.append((float(rightI[c] - leftI[c])) / (rightX - leftX))
            if renderM[0] == 6:
                dNCurr = []
                for c in range(3):
                    dNCurr.append((float(rightN[c] - leftN[c])) / (rightX - leftX))
        else:
            dZCurr = 0
            dICurr = [0,0,0]
            dNCurr = [0,0,0]

        #while this y level is not filled
        x = leftX
        while x < rightX + 1:
            #if coordinates are in the bounds of the canvas
            if x >= 0 and x < CanvasWidth - 1 and y >= 0 and y < CanvasHeight - 1:
                #if there isnt a "closer" pixel
                if z < zBuffer[math.floor(x)][math.floor(y)]:
                    #if edge pixel, use edge color
                    if renderM[0] == 2 and (x <= leftX + edgeThickness or x >= rightX - edgeThickness or y <= topY + edgeThickness or y >= botY - edgeThickness):
                        w.create_line(math.floor(x), math.floor(y), math.floor(x) + 1, math.floor(y), fill=color)
                    #if not edge pixel, use poly color
                    else:
                        #if not shading
                        if renderM[0] < 4:
                            w.create_line(math.floor(x), math.floor(y), math.floor(x) + 1, math.floor(y), fill=pColor)
                        #if flat shading
                        elif renderM[0] == 4:
                            w.create_line(math.floor(x), math.floor(y), math.floor(x) + 1, math.floor(y), fill=flatColor)
                        #if gouraud shading
                        elif renderM[0] == 5:
                            w.create_line(math.floor(x), math.floor(y), math.floor(x) + 1, math.floor(y), fill=triColorHexCode(intsy))
                            for c in range(3):
                                intsy[c] += dICurr[c]
                        #if phong shading
                        elif renderM[0] == 6:
                            w.create_line(math.floor(x), math.floor(y), math.floor(x) + 1, math.floor(y),
                                          fill=triColorHexCode(calcAmbDiffSpec(Ia, Kd, Ks, Ip, normsy, L, 1, V, specIndex)))
                            for c in range(3):
                                normsy[c] += dNCurr[c]
                    #set new zbuffer
                    zBuffer[math.floor(x)][math.floor(y)] = z
            z += dZCurr
            x += 1
        #set x and z boundaries for next y level
        eIX += edgeT[i].dX
        eJX += edgeT[j].dX
        eIZ += edgeT[i].dZ
        eJZ += edgeT[j].dZ
        #set intensity or normal boundaries for next y level if needed
        if renderM[0] == 5:
            for c in range(3):
                eII[c] += edgeT[i].dI[c]
                eJI[c] += edgeT[j].dI[c]
        if renderM[0] == 6:
            for c in range(3):
                eIN[c] += edgeT[i].dN[c]
                eJN[c] += edgeT[j].dN[c]
        #edge case dX checking, incase it becomes too big for some reason
        #checks for sign of dX/dZ; if positive, then checks if current edge x boundary is bigger than its able to
        #for negative, checks if smaller than able to
        if edgeT[i].dX > 0:
            if eIX > edgeT[i].xEnd:
                eIX = edgeT[i].xEnd
        if edgeT[i].dX < 0:
            if eIX < edgeT[i].xEnd:
                eIX = edgeT[i].xEnd
        if edgeT[j].dX > 0:
            if eJX > edgeT[j].xEnd:
                eJX = edgeT[j].xEnd
        if edgeT[j].dX < 0:
            if eJX < edgeT[j].xEnd:
                eJX = edgeT[j].xEnd

        #if bottom of edge is reached and polygon still not filled, switch to next edge
        if y >= edgeT[i].yEnd and y < botY:
            i = k
            eIX = edgeT[i].xStart
            eIZ = edgeT[i].zStart
            eIXEnd = edgeT[i].xEnd
            if renderM[0] == 5:
                for c in range(3):
                    eII[c] = edgeT[i].iStart[c]
            if renderM[0] == 6:
                for c in range(3):
                    eIN[c] = edgeT[i].nStart[c]
            k += 1
        if y >= edgeT[j].yEnd and y < botY:
            j = k
            eJX = edgeT[j].xStart
            eJZ = edgeT[j].zStart
            eJXEnd = edgeT[j].xEnd
            if renderM[0] == 5:
                for c in range(3):
                    eJI[c] = edgeT[j].iStart[c]
            if renderM[0] == 6:
                for c in range(3):
                    eJN[c] = edgeT[j].nStart[c]
            k += 1
        y += 1

#changes the render mode
def renderMode(event):
    #wireframe
    if event.keysym == "1":
        if renderM[0] != 1:
            renderM[0] = 1
            w.delete(ALL)
            drawObjects(objectList, objectColorList, zBuffer)
    #poly fill + edge
    if event.keysym == "2":
        if renderM[0] != 2:
            renderM[0] = 2
            w.delete(ALL)
            drawObjects(objectList, objectColorList, zBuffer)
    #poly fill
    if event.keysym == "3":
        if renderM[0] != 3:
            renderM[0] = 3
            w.delete(ALL)
            drawObjects(objectList, objectColorList, zBuffer)
    #flatshade
    if event.keysym == "4":
        if renderM[0] != 4:
            renderM[0] = 4
            w.delete(ALL)
            drawObjects(objectList, objectColorList, zBuffer)
    #gourad
    if event.keysym == "5":
        if renderM[0] != 5:
            renderM[0] = 5
            w.delete(ALL)
            drawObjects(objectList, objectColorList, zBuffer)
    #phong
    if event.keysym == "6":
        if renderM[0] != 6:
            renderM[0] = 6
            w.delete(ALL)
            drawObjects(objectList, objectColorList, zBuffer)

# gets edge between two points
# returns the two coordinates representing an edge
def getEdge(start,end,color,startNorm, endNorm):
    #project 3d to 2d, then convert projection to display coord
    startDisplay = convertToDisplayCoordinates(project(start))
    endDisplay = convertToDisplayCoordinates(project(end))
    #if shading get normal for point
    if renderM[0] > 4:
        startDisplay.append(startNorm)
        endDisplay.append(endNorm)
    return [startDisplay, endDisplay]

# Project the 3D endpoints to 2D point using a perspective projection implemented in 'project'
# Convert the projected endpoints to display coordinates via a call to 'convertToDisplayCoordinates'
# draw the actual line using the built-in create_line method
def drawLine(start,end,color):
    #project 3d to 2d, then convert projection to display coord
    startDisplay = convertToDisplayCoordinates(project(start))
    endDisplay = convertToDisplayCoordinates(project(end))
    w.create_line(startDisplay[0],startDisplay[1],endDisplay[0],endDisplay[1], fill=color)


# This function converts from 3D to 2D (+ depth) using the perspective projection technique.  Note that it
# will return a NEW list of points.  We will not want to keep around the projected points in our object as
# they are only used in rendering
def project(point):
    ps = []
    #conversion of 3d point to 2d point + depth
    for i in range(3):
        ps.append(d*(point[i]/(d+point[2])))
    return ps

# This function converts a 2D point to display coordinates in the tk system.  Note that it will return a
# NEW list of points.  We will not want to keep around the display coordinate points in our object as 
# they are only used in rendering.
def convertToDisplayCoordinates(point):
    displayCoord = []
    # conversion from 2d point to display coord
    displayCoord.append((CanvasWidth/2)+point[0])
    displayCoord.append((CanvasHeight/2)-point[1])
    displayCoord.append(point[2])
    return displayCoord
    
#gets reference point
# NEW list of points.  We will not want to keep around the display coordinate points in our object as
# they are only used in rendering.
def getRef(object):
    #get max and min for xyz
    refPoint = []
    negRefPoint = []
    minX = minY = minZ = sys.maxsize
    maxX = maxY = maxZ = -sys.maxsize - 1
    #for each vertex in the object(pyramid point cloud)
    for i in range(len(object)):
        #get max and min of each dimension
        minX = min(object[i][0], minX)
        minY = min(object[i][1], minY)
        minZ = min(object[i][2], minZ)
        maxX = max(object[i][0], maxX)
        maxY = max(object[i][1], maxY)
        maxZ = max(object[i][2], maxZ)
    #calc ref point
    refPoint.append((maxX+minX)/2)
    refPoint.append((maxY+minY)/2)
    refPoint.append((maxZ+minZ)/2)
    #negative stuff
    #dont know how python handles a negative 0 so i check for it
    #it may be unnecessary but eh
    for i in range(len(refPoint)):
        if refPoint[i] != 0:
            negRefPoint.append(-refPoint[i])
        else:
            negRefPoint.append(0)
    return refPoint, negRefPoint

#calculates 3d reflection vector
#takes a surface normal and a lighting vector
def reflect(N, L):
    R = []
    #normalize surface normal and lighting vector
    N = normalize(N)
    L = normalize(L)
    #le epic cosine trick
    twoCosPhi = 2  * (N[0]*L[0] + N[1]*L[1] + N[2]*L[2])
    #case checking
    if twoCosPhi > 0:
        for i in range(3):
            R.append(N[i] - (L[i] / twoCosPhi))
    elif twoCosPhi == 0:
        for i in range(3):
            R.append( -L[i])
    else:
        for i in range(3):
            R.append( -N[i] + (L[i] / twoCosPhi))
    return normalize(R)

#normalizes vectors
#math thing
#take a vector, spit a vector
def normalize(vector):
    sumOfSquares = 0
    for i in range(len(vector)):
        sumOfSquares += vector[i]**2
    magnitude = math.sqrt(sumOfSquares)
    vect = []
    for i in range(len(vector)):
        vect.append(vector[i]/magnitude)
    return vect

#calculates the ambient, diffuse, and specular lighting
def calcAmbDiffSpec(Ia, Kd, Ks, Ip, N, L, d, V, specIndex):
    #normalize normalize normalize
    L = normalize(L)
    V = normalize(V)
    #calc ambient lighting
    ambient = Ia * Kd
    #crossprod of surface normal and lighting vector
    NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2]
    if NdotL < 0:
        NdotL = 0
    #calc defuse
    diffuse = Ip * Kd * NdotL

    #reflection vector
    R = normalize(reflect(N,L))
    RdotV = R[0] * V[0] + R[1] * V[1] + R[2] * V[2]
    if RdotV < 0: RdotV = 0
    #calc specular
    specular = Ip * Ks * RdotV**specIndex
    return [ambient, diffuse, specular]

#turns colors into hexcodes
def colorHexCode(intensity):
    decInt = round(255 * intensity)
    if decInt > 255:
        print(decInt)
        decInt = 255
    hexString = str(hex(decInt))
    if hexString[0] == "-":
        trimmedHexString = "00"
    else:
        trimmedHexString = hexString[2:]
        if len(trimmedHexString) == 1:
            trimmedHexString = "0" + trimmedHexString
    return trimmedHexString

#combines 3 hex color codes into one whole rgb code
def triColorHexCode(ambDiffSpec):
    combinedColorCode = colorHexCode(ambDiffSpec[0]+ambDiffSpec[1]+ambDiffSpec[2])
    specularColorCode = colorHexCode(ambDiffSpec[2])
    colorString = "#" + specularColorCode + combinedColorCode + specularColorCode
    return colorString


# **************************************************************************
# Everything below this point implements the interface
def reset():
    w.delete(ALL)
    resetPyramid()
    drawObjects(objectList, objectColorList, zBuffer)

def larger():
    w.delete(ALL)
    scale(objectPointCloudList[currentObj[0]],1.1)
    drawObjects(objectList, objectColorList, zBuffer)

def smaller():
    w.delete(ALL)
    scale(objectPointCloudList[currentObj[0]],.9)
    drawObjects(objectList, objectColorList, zBuffer)

def forward():
    w.delete(ALL)
    translate(objectPointCloudList[currentObj[0]],[0,0,5])
    drawObjects(objectList, objectColorList, zBuffer)

def backward():
    w.delete(ALL)
    translate(objectPointCloudList[currentObj[0]],[0,0,-5])
    drawObjects(objectList, objectColorList, zBuffer)

def left():
    w.delete(ALL)
    translate(objectPointCloudList[currentObj[0]],[-5,0,0])
    drawObjects(objectList, objectColorList, zBuffer)

def right():
    w.delete(ALL)
    translate(objectPointCloudList[currentObj[0]],[5,0,0])
    drawObjects(objectList, objectColorList, zBuffer)

def up():
    w.delete(ALL)
    translate(objectPointCloudList[currentObj[0]],[0,5,0])
    drawObjects(objectList, objectColorList, zBuffer)

def down():
    w.delete(ALL)
    translate(objectPointCloudList[currentObj[0]],[0,-5,0])
    drawObjects(objectList, objectColorList, zBuffer)

def xPlus():
    w.delete(ALL)
    rotateX(objectPointCloudList[currentObj[0]],5)
    drawObjects(objectList, objectColorList, zBuffer)

def xMinus():
    w.delete(ALL)
    rotateX(objectPointCloudList[currentObj[0]],-5)
    drawObjects(objectList, objectColorList, zBuffer)

def yPlus():
    w.delete(ALL)
    rotateY(objectPointCloudList[currentObj[0]],5)
    drawObjects(objectList, objectColorList, zBuffer)

def yMinus():
    w.delete(ALL)
    rotateY(objectPointCloudList[currentObj[0]],-5)
    drawObjects(objectList, objectColorList, zBuffer)

def zPlus():
    w.delete(ALL)
    rotateZ(objectPointCloudList[currentObj[0]],5)
    drawObjects(objectList, objectColorList, zBuffer)

def zMinus():
    w.delete(ALL)
    rotateZ(objectPointCloudList[currentObj[0]],-5)
    drawObjects(objectList, objectColorList, zBuffer)


def prevObj():
    w.delete(ALL)
    changeObj(-1)
    drawObjects(objectList, objectColorList, zBuffer)

def nextObj():
    w.delete(ALL)
    changeObj(1)
    drawObjects(objectList, objectColorList, zBuffer)

root = Tk()
outerframe = Frame(root)
outerframe.pack()

w = Canvas(outerframe, width=CanvasWidth, height=CanvasHeight)
drawObjects(objectList, objectColorList, zBuffer)
w.pack()

controlpanel = Frame(outerframe)
controlpanel.pack()

resetcontrols = Frame(controlpanel, height=100, borderwidth=2, relief=RIDGE)
resetcontrols.pack(side=LEFT)

resetcontrolslabel = Label(resetcontrols, text="Reset")
resetcontrolslabel.pack()

resetButton = Button(resetcontrols, text="Reset", fg="green", command=reset)
resetButton.pack(side=LEFT)

scalecontrols = Frame(controlpanel, borderwidth=2, relief=RIDGE)
scalecontrols.pack(side=LEFT)

scalecontrolslabel = Label(scalecontrols, text="Scale")
scalecontrolslabel.pack()

largerButton = Button(scalecontrols, text="Larger", command=larger)
largerButton.pack(side=LEFT)

smallerButton = Button(scalecontrols, text="Smaller", command=smaller)
smallerButton.pack(side=LEFT)

translatecontrols = Frame(controlpanel, borderwidth=2, relief=RIDGE)
translatecontrols.pack(side=LEFT)

translatecontrolslabel = Label(translatecontrols, text="Translation")
translatecontrolslabel.pack()

forwardButton = Button(translatecontrols, text="FW", command=forward)
forwardButton.pack(side=LEFT)

backwardButton = Button(translatecontrols, text="BK", command=backward)
backwardButton.pack(side=LEFT)

leftButton = Button(translatecontrols, text="LF", command=left)
leftButton.pack(side=LEFT)

rightButton = Button(translatecontrols, text="RT", command=right)
rightButton.pack(side=LEFT)

upButton = Button(translatecontrols, text="UP", command=up)
upButton.pack(side=LEFT)

downButton = Button(translatecontrols, text="DN", command=down)
downButton.pack(side=LEFT)

rotationcontrols = Frame(controlpanel, borderwidth=2, relief=RIDGE)
rotationcontrols.pack(side=LEFT)

rotationcontrolslabel = Label(rotationcontrols, text="Rotation")
rotationcontrolslabel.pack()

xPlusButton = Button(rotationcontrols, text="X+", command=xPlus)
xPlusButton.pack(side=LEFT)

xMinusButton = Button(rotationcontrols, text="X-", command=xMinus)
xMinusButton.pack(side=LEFT)

yPlusButton = Button(rotationcontrols, text="Y+", command=yPlus)
yPlusButton.pack(side=LEFT)

yMinusButton = Button(rotationcontrols, text="Y-", command=yMinus)
yMinusButton.pack(side=LEFT)

zPlusButton = Button(rotationcontrols, text="Z+", command=zPlus)
zPlusButton.pack(side=LEFT)

zMinusButton = Button(rotationcontrols, text="Z-", command=zMinus)
zMinusButton.pack(side=LEFT)

#i hate tkinter
#buttons for object selection
objectSelection = Frame(controlpanel, borderwidth=2, relief=RIDGE)
objectSelection.pack(side=LEFT)

objectSelectionlabel = Label(objectSelection, text="Object Selection")
objectSelectionlabel.pack()

prevObj = Button(objectSelection, text="Previous", command=prevObj)
prevObj.pack(side=LEFT)
nextObj = Button(objectSelection, text="Next", command=nextObj)
nextObj.pack(side=LEFT)

root.bind("<Key>", renderMode)

root.mainloop()
