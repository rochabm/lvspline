# Parametrized Left Ventricle Geometry and Mesh Generator
# Gmsh is used as the backend
#
# Authors
# Bernardo Martins Rocha <bernardomartinsrocha@gmail.com>
# Joventino de Oliveira Campos <joventinoo@gmail.com>
# January, 2019

import sys, os
import numpy as np
from math import sin, cos, pi, tan, asin, acos, degrees

# Variaveis globais
malhaTet = True

# GMSH STUFF
gmsh_line_id = 0
gmsh_vertex_id = 0
gmsh_surface_id = 0
gmsh_volume_id = 0
refinement_id = 1

geoFile = None

#Angle of the curvature in the base of the ellpsoid. (original value =  pi/6)
curvature_angle = 0.00001 #pi/6

# ------------------------------------------------------------------------------

class Block:
    """
    Holds data from one segment (block) from the AHA-17 segments model of the LV
    """
    def __init__(self, id, a1_1, a1_2, a1_3, a1_4, a2_1, a2_2, a2_3, a2_4, t0, t1, p0, p1):
        self.num = id
        self.phi0 = p0
        self.phi1 = p1
        self.theta0 = t0
        self.theta1 = t1
        self.apex = False

        if (a1_3 == None and a1_4 == None and a2_3 == None and a2_4 == None):
            a1_3, a1_4, a2_3, a2_4 = -1, -1, -1, -1
            self.apex = True

        self.a_endo1, self.a_endo2 = a1_1, a1_2
        self.a_endo3, self.a_endo4 = a1_3, a1_4
        self.a_epi1, self.a_epi2 = a2_1, a2_2
        self.a_epi3, self.a_epi4 = a2_3, a2_4

    def __str__(self):
        a, b, c, d = self.a_endo1, self.a_endo2, self.a_endo3, self.a_endo4
        e, f, g, h = self.a_epi1, self.a_epi2, self.a_epi3, self.a_epi4
        s = 'Block %d ;\n' % self.num
        s += '  a_endo = (%f, %f, %f, %f);\n' % (a, b, c, d)
        s += '  a_epi  = (%f, %f, %f, %f);\n' % (e, f, g, h)
        s += '  phi    = (%f, %f);\n' % (self.phi0, self.phi1)
        s += '  theta  = (%f, %f);\n' % (self.theta0, self.theta1)
        s += '  apex   = %s' % self.apex
        return s

    def getNum(self):
        return int(self.num)

    def getEndoData(self):
        a11, a12, a13, a14 = self.a_endo1, self.a_endo2, self.a_endo3, self.a_endo4
        return a11, a12, a13, a14

    def getEpiData(self):
        a21, a22, a23, a24 = self.a_epi1, self.a_epi2, self.a_epi3, self.a_epi4
        return a21, a22, a23, a24

    def getTheta(self):
        return self.theta0, self.theta1

    def getPhi(self):
        return self.phi0, self.phi1

    def inside(self, teta, phi):
        if (phi > 0):
            if (teta >= self.theta0 and teta <= self.theta1 and
                        phi >= self.phi0 and phi <= self.phi1):
                return True
        elif (phi < 0):
            if (teta >= self.theta0 and teta <= self.theta1 and
                        phi >= self.phi0 and phi <= self.phi1):
                return True
        else:
            return False

    def isApex(self):
        return self.apex

# ------------------------------------------------------------------------------

def findBlock(blocklist, teta, phi):
    myblock, blockid = None, None
    for i in range(len(blocklist)):
        res = blocklist[i].inside(teta, phi)
        if (res):
            myblock = blocklist[i]
            blockid = i + 1
            #print i
    if (blockid is None):
        print("ERROR BLOCK NOT FOUND")
        
        return None
    return myblock

# ------------------------------------------------------------------------------

def setCharacteristicLength(value):
    global geoFile
    geoFile.write("lc = %f;\n" % (value))

#------------------------------------------------------------------------------

def createVertex(x, y, z):
    global gmsh_vertex_id, geoFile
    gmsh_vertex_id += 1
    geoFile.write("Point(%d) = {%e, %e, %e, lc};\n" % (gmsh_vertex_id,x,y,z))
    return gmsh_vertex_id

# ------------------------------------------------------------------------------

def createLine(nodes):
    global gmsh_line_id, geoFile
    gmsh_line_id += 1
    geoFile.write("Line(%d) = {%d, %d};\n" % (gmsh_line_id,nodes[0],nodes[1]))
    return gmsh_line_id

# ------------------------------------------------------------------------------

def createSpline(vertexList):
    global gmsh_line_id, geoFile
    gmsh_line_id += 1
    cmdstr = 'Spline(%d) = {' % gmsh_line_id
    for i in range(len(vertexList)-1):
        cmdstr += ' %d,' % vertexList[i]
    cmdstr += ' %d' % vertexList[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)
    return gmsh_line_id

# ------------------------------------------------------------------------------

def createCompoundSurface(surfID):
    global gmsh_surface_id, geoFile
    
    gmsh_surface_id +=1

    cmdstr = 'Compound Surface(%d) = {' % gmsh_surface_id
    for i in range(len(surfID)-1):
        cmdstr += ' %d,' % surfID[i]
    cmdstr += ' %d' % surfID[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)

    return gmsh_surface_id

# ------------------------------------------------------------------------------

def createPhysicalSurface(surfID, marker):
    global geoFile

    if isinstance(surfID,int):
        geoFile.write('Physical Surface(%d) = {%d};\n' % (marker, surfID) )
    else:
        cmdstr = 'Physical Surface(%d) = {' % marker
        for i in range(len(surfID)-1):
            cmdstr += ' %d,' % surfID[i]
        cmdstr += ' %d' % surfID[-1]
        cmdstr += ' };\n'
        geoFile.write("%s" % cmdstr)

# ------------------------------------------------------------------------------

def createPhysicalVolume(volID, marker):
    global geoFile
    geoFile.write('Physical Volume(%d) = {%d};\n' % (marker, volID) )

# ------------------------------------------------------------------------------

def createSurfaceBase(endoBase, epiBase):
    global gmsh_line_id, gmsh_surface_id, geoFile
  
    gmsh_surface_id += 1

    cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
    for i in range(len(endoBase)-1):
        cmdstr += ' %d,' % endoBase[i]

    cmdstr += ' %d' % endoBase[-1] 
    cmdstr += ' };\n'

    geoFile.write("%s" % cmdstr)

    cmdstr2 = "Plane Surface(%d) = {%d" % (gmsh_surface_id,gmsh_surface_id)
    bID = gmsh_surface_id
    
    gmsh_surface_id += 1
    cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
    for i in range(len(epiBase)-1):
        cmdstr += ' %d,' % epiBase[i]
    
    cmdstr += ' %d' % epiBase[-1]  
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)

    cmdstr2+=", %d};\n" % (gmsh_surface_id)
    geoFile.write("%s" % cmdstr2)

    gmsh_surface_id += 1
    geoFile.write('Surface Loop(%d) = {%d};\n' % (gmsh_surface_id, bID) )

    #Marking epicardium line
    cmdstr = 'Physical Line(50) = {'
    for i in range(len(epiBase)-1):
        cmdstr += ' %d,' % epiBase[i]

    cmdstr += ' %d' % epiBase[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)

    return gmsh_surface_id, bID

# ------------------------------------------------------------------------------

def createSurfaceSkin(splineList, splineVertex):

    global gmsh_line_id, gmsh_surface_id, geoFile

    baseSplines1 = []
    baseSplines2 = []
    baseSplines3 = []
    baseSplines4 = []
    surfaceLoop = []

    for i in range(len(splineList[0])-1):

        sid1 = createSpline([splineVertex[0][i], splineVertex[0][i+1]])
        baseSplines1.append(sid1)
        sid2 = createSpline([splineVertex[1][i], splineVertex[1][i+1]])
        baseSplines2.append(sid2)
        sid3 = createSpline([splineVertex[2][i], splineVertex[2][i+1]])
        baseSplines3.append(sid3)
        sid4 = createSpline([splineVertex[3][i], splineVertex[3][i+1]])
        baseSplines4.append(sid4)
        sid = [sid1, sid2, sid3, sid4]

        for j in range(3):
            gmsh_surface_id += 1
            cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
            cmdstr += ' %d,' % sid[j]
            cmdstr += ' %d,' % splineList[j][i+1]
            cmdstr += ' %d,' % -sid[j+1]
            cmdstr += ' %d' % -splineList[j][i]
            cmdstr += ' };\n'
            geoFile.write("%s" % cmdstr)
            geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
            surfaceLoop.append(gmsh_surface_id)


        gmsh_surface_id += 1
        cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
        cmdstr += ' %d,' % sid[3]
        cmdstr += ' %d,' % splineList[3][i+1]
        cmdstr += ' %d' % -splineList[3][i]
        cmdstr += ' };\n'
        geoFile.write("%s" % cmdstr)
        geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
        surfaceLoop.append(gmsh_surface_id)

    sid1 = createSpline([splineVertex[0][-1], splineVertex[0][0]])
    baseSplines1.append(sid1)
    sid2 = createSpline([splineVertex[1][-1], splineVertex[1][0]])
    baseSplines2.append(sid2)
    sid3 = createSpline([splineVertex[2][-1], splineVertex[2][0]])
    baseSplines3.append(sid3)
    sid4 = createSpline([splineVertex[3][-1], splineVertex[3][0]])
    baseSplines4.append(sid4)
    sid = [sid1, sid2, sid3, sid4]

    for j in range(3):
        gmsh_surface_id += 1
        cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
        cmdstr += ' %d,' % sid[j]
        cmdstr += ' %d,' % splineList[j][0]
        cmdstr += ' %d,' % -sid[j+1]
        cmdstr += ' %d' % -splineList[j][-1]
        cmdstr += ' };\n'
        geoFile.write("%s" % cmdstr)
        geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
        surfaceLoop.append(gmsh_surface_id)

    gmsh_surface_id += 1
    cmdstr = 'Line Loop(%d) = {' % gmsh_surface_id
    cmdstr += ' %d,' % sid[3]
    cmdstr += ' %d,' % splineList[3][0]
    cmdstr += ' %d' % -splineList[3][-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)
    geoFile.write("Ruled Surface(%d) = {%d};\n" % (gmsh_surface_id,gmsh_surface_id))
    surfaceLoop.append(gmsh_surface_id)

    surfID = gmsh_surface_id

    gmsh_surface_id +=1
    cmdstr = 'Surface Loop(%d) = {' % gmsh_surface_id
    for i in range(len(surfaceLoop)-1):
        cmdstr += ' %d,' % surfaceLoop[i]
    cmdstr += ' %d' % surfaceLoop[-1]
    cmdstr += ' };\n'
    geoFile.write("%s" % cmdstr)

    return gmsh_surface_id, surfaceLoop, baseSplines1

# ------------------------------------------------------------------------------

def createVolume(surfList):
    endoSID, epiSID, baseSID = surfList
    global gmsh_volume_id, gmsh_surface_id, geoFile
    gmsh_volume_id +=1

    cmdstr = 'Volume(%d) = {%d, %d, %d};\n' % (gmsh_volume_id, endoSID, baseSID, epiSID)
    geoFile.write("%s" % cmdstr)

    return gmsh_volume_id

# ------------------------------------------------------------------------------

def ellipsoid(a, c, phi, theta):
    """
    Parametric equations for the ellipsoid
    """
    #translate base to 0
    offset = c * sin(curvature_angle)
    x = a * cos(phi) * cos(theta)
    y = a * cos(phi) * sin(theta)
    z = c * sin(phi) + offset
    return (x, y, z)

# ------------------------------------------------------------------------------

def lerp(r0, r1, theta0, theta1, theta):
    """
    Linear interpolation
    """
    r = r0 + ((r1 - r0) / (theta1 - theta0)) * (theta - theta0)
    return r

# -------------------------------------------------------------------------------

def createBlocks(a1, c1, c2, bull):
    nblocks = 0
    lblocks = []

    # increments according the ellipsoid parametrization
    # z = c * sin (phi)
    # phi = asin(z/c)
    inc1 = asin(1./3.) # basal region angle
    inc2 = asin(2./3.) # medial region angle
    inc3 = asin(3./3.) # apical region angle

    # AQUI
    phi0 = -curvature_angle
    phi1 = -curvature_angle + inc1
    phi2 = -curvature_angle + inc2
    phi3 = -curvature_angle + inc3

    print("Generating blocks in basal region")
    for i in range(6):
        nblocks += 1
        print("Block %d" % (nblocks))
        theta0 = (0 + i) * pi / 3.0
        theta1 = (1 + i) * pi / 3.0
        if (i < 5):
            a2_0 = a1 + (bull[0 + i] );
            a2_1 = a1 + (bull[1 + i] )
            a2_2 = a1 + (bull[6 + i] );
            a2_3 = a1 + (bull[7 + i] )
            a1_0 = a1 
            a1_1 = a1 
            a1_2 = a1 
            a1_3 = a1 

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi0, phi1)
            lblocks.append(b)
        else:
            a2_0 = a1 + (bull[5 + (i - 5)] );
            a2_1 = a1 + (bull[0 + (i - 5)] )
            a2_2 = a1 + (bull[11 + (i - 5)] );
            a2_3 = a1 + (bull[6 + (i - 5)] )
            a1_0 = a1 
            a1_1 = a1 
            a1_2 = a1 
            a1_3 = a1 

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi0, phi1)
            lblocks.append(b)

    print("Generating blocks in medial region")
    for i in range(6):
        nblocks += 1
        print("Block %d" % (nblocks))
        theta0 = (0 + i) * pi / 3
        theta1 = (1 + i) * pi / 3
        if (i < 5):
            a2_0 = a1 + (bull[i + 6] );
            a2_1 = a1 + (bull[i + 7] )
            if (i==0):
                raio0 = lerp(a1 + bull[12], a1 + bull[13], 0, pi/2, pi/12)
                a2_2 = raio0
                raio1 = lerp(a1 + bull[12], a1 + bull[13], 0, pi/2, pi/2-pi/12)
                a2_3 = raio1
            elif (i==1):
                raio0 = lerp(a1 + bull[12], a1 + bull[13], 0, pi/2, pi/2-pi/12)
                a2_2 = raio0
                raio1 = lerp(a1 + bull[13], a1 + bull[14], 0, pi/2, 2*pi/3 - (pi/2-pi/12))
                a2_3 = raio1
            elif (i==2):
                raio0 = lerp(a1 + bull[13], a1 + bull[14], 0, pi/2, 2*pi/3 - (pi/2-pi/12))
                a2_2 = raio0
                raio1 = lerp(a1 + bull[14], a1 + bull[15], 0, pi/2, pi/12)
                a2_3 = raio1
            elif (i==3):
                raio0 = lerp(a1 + bull[14], a1 + bull[15], 0, pi/2, pi/12)
                a2_2 = raio0
                raio1 = lerp(a1 + bull[14], a1 + bull[15], 0, pi/2, 4*pi/3 - (pi - pi/12) )
                a2_3 = raio1
            elif (i==4):
                raio0 = lerp(a1 + bull[14], a1 + bull[15], 0, pi/2, 4*pi/3 - (pi - pi/12) )
                a2_2 = raio0
                raio1 = lerp(a1 + bull[15], a1 + bull[16], 0, pi/2, 5*pi/3 - (3*pi/2 - pi/12) )
                a2_3 = raio1


            a1_0 = a1 
            a1_1 = a1 
            a1_2 = a1 
            a1_3 = a1 

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi1, phi2)
            lblocks.append(b)
        else:
            a2_0 = a1 + (bull[5 + (i + 6 - 5)] )
            a2_1 = a1 + (bull[0 + (i + 6 - 5)] )

            raio0 = lerp(a1 + bull[15], a1 + bull[16], 0, pi/2, 5*pi/3 - (3*pi/2 - pi/12) )
            a2_2 = raio0 
            raio1 = lerp(a1 + bull[12], a1 + bull[13], 0, pi/2, pi/12)
            a2_3 = raio1

            a1_0 = a1 
            a1_1 = a1 
            a1_2 = a1 
            a1_3 = a1 

            b = Block(nblocks, a1_0, a1_1, a1_2, a1_3, a2_0, a2_1, a2_2, a2_3,
                      theta0, theta1, phi1, phi2)
            lblocks.append(b)

    
    print("Generating blocks in apical region")
    for i in range(5):
        nblocks += 1
        print("Block %d" % (nblocks))
        block = i
        raio0 = lerp(a1 + bull[12], a1 + bull[13], 0, pi/2, pi/12)
        #-pi/12 is the rotation in apical section of the Bulls eye diagram
        if (i==0):
            theta0 = 0
        else:
            theta0 = (0 + i) * pi / 2 - pi/12
        theta1 = (1 + i) * pi / 2 - pi/12
        if (i < 4):
            a2_0 = a1 + (bull[i + 12] );
            a2_1 = a1 + (bull[i + 13] )
            a1_0 = a1 
            a1_1 = a1 
            
            b = Block(nblocks, a1_0, a1_1, None, None, a2_0, a2_1, None, None,
                      theta0, theta1, phi2, phi3)
            lblocks.append(b)
        else:
            a2_0 = a1 + (bull[16] );
            a2_1 = raio0 
            a1_0 = a1 
            a1_1 = a1 

            b = Block(nblocks, a1_0, a1_1, None, None, a2_0, a2_1, None, None,
                      -pi/12 + 2*pi, 2*pi, phi2, phi3)
            lblocks.append(b)
    
        

    print("Blocks generation OK")
    return lblocks

# ------------------------------------------------------------------------------

def getRegion(x, y, z, a, c, blocklist):
    offset = c * sin(curvature_angle)
    teste_phi = asin((z - offset)/c)
    teste_theta = asin( y/(a*cos(teste_phi)) )
    if teste_theta < 0: 
        teste_theta = teste_theta + 2*pi

    myblock = findBlock(blocklist, teste_theta, teste_phi)

    return myblock.getNum()

# ------------------------------------------------------------------------------

def geraLVSplines(bull, blocks, a1, c1, c2, phi0):

    # initilization
    ncirc = 18*2
    nlong = 18

    phi0, phi1 = -curvature_angle, pi / 2
    theta0, theta1 = 0, 2 * pi

    dt = (theta1) / (ncirc - 1)
    dp = (phi1 - phi0) / (nlong - 1)

    
    print("Gmsh .geo file generation;\n;\n;\n")

    # create apex nodes - endo and epi
    x, y, z = ellipsoid(0, c1, phi1, 0)
    apexEndo = createVertex(x, y, z)

    x, y, z = ellipsoid(0, c2, phi1, 0)
    apexEpi = createVertex(x, y, z)

    LEpiBase = []
    LEndoBase = []
    ZEndoBase = []

    #
    # epicardium
    #
    theta = theta0
    
    splst1 = []
    splst2= []
    splst3 = []
    splst4 = []
    EpisplstVertex1 = []
    EpisplstVertex2 = []
    EpisplstVertex3 = []
    EpisplstVertex4 = []

    for j in range(ncirc-1):
        phi = phi0
        vlist = []
        for i in range(nlong - 1):
            myblock = findBlock(blocks, theta, phi)
            
            a1_0, a1_1, a1_2, a1_3 = myblock.getEpiData()

            bphi0, bphi1 = myblock.getPhi()
            btet0, btet1 = myblock.getTheta()

            if (not myblock.isApex()):
                ai = lerp(a1_0, a1_2, bphi0, bphi1, phi)
                af = lerp(a1_1, a1_3, bphi0, bphi1, phi)
                a = lerp(ai, af, btet0, btet1, theta)
            else:
                a = lerp(a1_0, a1_1, btet0, btet1, theta)

            x, y, z = ellipsoid(a, c2, phi, theta)

            
            vid = createVertex(x, y, z)


            vlist.append(vid)

            if phi == phi0:
                LEpiBase.append(vid)
                ZEndoBase.append(z)

            phi = phi + dp

        vlist.append(apexEpi)

        sid = createSpline(vlist[0:7])
        splst1.append(sid)
        EpisplstVertex1.append(vlist[0])
        sid = createSpline(vlist[6:11])
        splst2.append(sid)
        EpisplstVertex2.append(vlist[6])
        sid = createSpline(vlist[10:15])
        splst3.append(sid)
        EpisplstVertex3.append(vlist[10])
        sid = createSpline(vlist[14:])
        splst4.append(sid)
        EpisplstVertex4.append(vlist[14])


        theta = theta + dt

    # create epicardium surface
    epiSID, surfID, epiBaseSplines = createSurfaceSkin([splst1, splst2, splst3, splst4],
                                                       [EpisplstVertex1, EpisplstVertex2,
                                                        EpisplstVertex3, EpisplstVertex4])
    ID=createCompoundSurface(surfID)
    createPhysicalSurface(ID, 40)
    
    #
    # endocardium 
    #
    theta = theta0
    splst1 = []
    splst2= []
    splst3 = []
    splst4 = []
    EndosplstVertex1 = []
    EndosplstVertex2 = []
    EndosplstVertex3 = []
    EndosplstVertex4 = []

    for j in range(ncirc - 1):
        phi = phi0
        vlist = []
        for i in range(nlong - 1):
            # find current block
            myblock = findBlock(blocks, theta, phi)

            # get block information
            a1_0, a1_1, a1_2, a1_3 = myblock.getEndoData()

            bphi0, bphi1 = myblock.getPhi()
            btet0, btet1 = myblock.getTheta()

            if (not myblock.isApex()):
                ai = lerp(a1_0, a1_2, bphi0, bphi1, phi)
                af = lerp(a1_1, a1_3, bphi0, bphi1, phi)
                a = lerp(ai, af, btet0, btet1, theta)
            else:
                a = lerp(a1_0, a1_1, btet0, btet1, theta)

            # create point
            x, y, z = ellipsoid(a, c1, phi, theta)        
            
            if phi == phi0: z=ZEndoBase[j]

            
            vid = createVertex(x, y, z)

            # add vertex to the list
            vlist.append(vid)

            if phi == phi0: 
                LEndoBase.append(vid)

            # increment phi angle
            phi = phi + dp

        vlist.append(apexEndo)

        
        # add spline in the list
        sid = createSpline(vlist[0:7])
        splst1.append(sid)
        EndosplstVertex1.append(vlist[0])
        sid = createSpline(vlist[6:11])
        splst2.append(sid)
        EndosplstVertex2.append(vlist[6])
        sid = createSpline(vlist[10:15])
        splst3.append(sid)
        EndosplstVertex3.append(vlist[10])
        sid = createSpline(vlist[14:])
        splst4.append(sid)
        EndosplstVertex4.append(vlist[14])
        
        # increment theta angle
        theta = theta + dt
        
    # create endocardium surface
    endoSID, surfID, endoBaseSplines = createSurfaceSkin([splst1, splst2, splst3, splst4],
                                                         [EndosplstVertex1, EndosplstVertex2,
                                                          EndosplstVertex3, EndosplstVertex4])
    ID = createCompoundSurface(surfID)
    createPhysicalSurface(ID, 30)

    # defines splines at the base
    baseSID, surfID = createSurfaceBase(endoBaseSplines, epiBaseSplines)
    createPhysicalSurface(surfID, 10)

    # surfaces ids
    surflist = [endoSID, epiSID, baseSID]
    veID = createVolume(surflist)
    createPhysicalVolume(veID, 0)

# -----------------------------------------------------------------------------

def createTetMesh(baseName):
    cmdstr = 'gmsh ' + baseName + '.geo -3 -o ' + baseName + '.msh'
    os.system(cmdstr)

# ------------------------------------------------------------------------------

def refineMesh():
    global geoFile, refinement_id
    refinement_id +=1
    geoFile.write('Mesh %d;\n' % (refinement_id))
    geoFile.write('RefineMesh;\n')
    
# -----------------------------------------------------------------------------

def createFibers(baseName, params):
    
    #from dolfin import *
    #from fiberrules import *
    #from dolfin_utils.meshconvert import meshconvert
    
    parameters["form_compiler"]["quadrature_degree"] = 2
    parameters.allow_extrapolation = True
    
    #convert mesh to dolfin format in order to use in Fiberrules
    ifilename = baseName + '.msh'
    ofilename = baseName + '_fenics.xml'
    iformat = 'gmsh'
    meshconvert.convert2xml(ifilename, ofilename, iformat=iformat)
    
    #Load mesh for fiberrules
    mesh = Mesh(ofilename)
    domains = MeshFunction("size_t", mesh, baseName + '_fenics_facet_region.xml')
    
    for facet in facets(mesh):
        mesh.domains().set_marker((facet.index(), domains[facet]), 2)
       
    fiber_angle_epi = params['fiber_angle_epi']
    fiber_angle_endo = params['fiber_angle_endo']
    sheet_angle_endo = params['sheet_angle_endo']
    sheet_angle_epi = params['sheet_angle_epi']
    
    fiber_space = FunctionSpace(mesh, 'DG',0)
    
    fib, she, crs = dolfin_fiberrules(mesh, fiber_space, fiber_angle_epi, fiber_angle_endo,
                                      sheet_angle_epi, sheet_angle_endo)

    dolfin_to_vtk(fib, baseName)
    dolfin_to_carpfile(mesh, baseName)
    if params['fiber_type'] == "transversely_isotropic":
        fibers_to_carpfile(mesh, fib, baseName)
    elif params['fiber_type'] == "orthotropic":
        fibers_to_carpfile_ortho(mesh, fib, she, crs, baseName)
    else:
        print('Error: fiber type not defined!')
        exit(0)

# ------------------------------------------------------------------------------

def generateLVMesh(a1, c1, original_bull, baseName='out_lv_mesh'):
    global geoFile

    # Some initial configs
    bull_ids = (0, 5, 4, 3, 2, 1, 6, 11, 10, 9, 8, 7, 12, 15, 14, 13, 12, 16)
    bull = []
    for i in bull_ids:
        bull.append(original_bull[i])    

    fileName = baseName + '.geo'
    geoFile = open(fileName,'w')

    c2 = c1 + bull[-1]
    
    # change long axis in order to consider the LV curvature 
    reduction = c1*sin(curvature_angle)/(1.0 + sin(curvature_angle))
    c1 = c1 - reduction
    a1 = a1/2.0 

    reduction = c2*sin(curvature_angle)/(1.0 + sin(curvature_angle))
    c2 = c2 - reduction

    # create list with the blocks and their information
    blocklist = createBlocks(a1, c1, c2, bull)

    # Define geometry
    setCharacteristicLength(1.0)
    geraLVSplines(bull, blocklist, a1, c1, c2, 0)
    #refineMesh()   
    geoFile.close()

    # Generates finite element mesh
    createTetMesh(baseName)

    # This code is not available (it depends on DOLFIN/fiberrules package)
    #createFibers(baseName)
    #convert2cardiax(...)

# End
