#Pygame is used for the main interface
import pygame
from pygame.locals import *

#Used as a bridge between pygame and opengl for 3D visualization
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np

import math

import trackingV5

#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!README:  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#   This is a visualizing tool to simulate a heliostat system that can move around
#   in a 3D plane.
#
#   Matrix/Vector nomenclature:
#   Local Axis:  Are the blue (Z), green (X), and red (Y) axis's that define
#   the heliostat itself.  These move with user inputted rotations (see below)
#
#   Absolute Visual Vectex:  The white axis, that are absolute relative to the user
#
#   Heliostat vertex:  These are to the heliostat (rectangle) itself.
#
#   Input Buttons:
#   q,w: Rotate around the local X axis
#   a,s: Rotate around the local Y axis
#   z,x: Rotate around the local Z axis
#   e: Rotate -90 around the local X axis:  I recommend using this as soon as the
#       screen starts to orient the Z in the up direction (towards the sky)
#   -,=: Zoom out, in
#   g,h: Turns on/off the ground (green) and the sky (blue)
#   r: turns on/off sun rays
#   p: Turn on/off the tracking algorithm

#Parameters for distances between power tower, heliostat
stationaryMirrorHeight = 12.0
stationaryMirrorWidth = 9.0
distanceStationaryMirrorToRotatedMirror = 32

#Dimensions of the heliostat mirror, and distances between mirrors
heliostatMirrorHeight = 12.0
heliostatMirrorWidth = 9.0
helioOneOffset = -18
helioTwoOffset = -7
helioThreeOffset = 7
helioFourOffset = 18
helioEnclosureX1, helioEnclosureX2 = 35.0, -21.0
helioEnclosureY = 8.0

powerTowerDimension = 12
powerTowerYoffset = -112
powerTowerZoffset = 63

#Output:  Draws the x,y,z axis on the screen
#Green = X axis
#Red = Y axis
#Blue = Z axis (Toward the sky)
def Axis(absYaxis,absXaxis):
    axis = (
        (0, 0, 0),
        (3, 0, 0),
        (0, 3, 0),
        (0, 0, 3),
        #These are the visual absolute coordinates
        (absYaxis[0],absYaxis[1],absYaxis[2]),
        (absXaxis[0],absXaxis[1],absXaxis[2])
    )
    axisEdge = (
        (0,1),
        (0,2),
        (0,3),
        (0,4),
        (0,5),
    )
    #Green (X), Red (Y), Blue (Z), white (2)
    axisColors = ((0,255,0),(255,0,0),(0,0,255),(255,255,255),(255,255,255))
    x = 0
    glBegin(GL_LINES)
    for edge in axisEdge:
        glColor3fv(axisColors[x])
        x = x + 1
        for vertex in edge:
            glVertex3fv(axis[vertex])
    glEnd()

#Input: degrees to rotate heliostat around power tower
#Output: Change in x and y distance to move the heliostat radially around the power tower
def moveHelio(degrees):
    radians = degrees * math.pi / 180
    y = -(powerTowerYoffset * math.cos(radians) - powerTowerYoffset)
    x = -powerTowerYoffset * math.sin(radians)
    return x, y

class rotatedMirror:
    def __init__(self, degrees):
        self.degrees = degrees
        #Used for mapping the coordinates of the corners of the moving heliostat

        #These are the local axis coordinates for all rotated mirrors
        #And the enclosure around the mirrors after the radial rotation
        #around the power tower
        self.helioCoord = np.array([
        #Heliostat Two
        np.array([
            np.array([heliostatMirrorWidth/2+helioTwoOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([heliostatMirrorWidth/2+helioTwoOffset, heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioTwoOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioTwoOffset, heliostatMirrorHeight/2, 0.0, 1]),
        ]),
        #Heliostat One
        np.array([
            np.array([heliostatMirrorWidth/2+helioOneOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([heliostatMirrorWidth/2+helioOneOffset, heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioOneOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioOneOffset, heliostatMirrorHeight/2, 0.0, 1]),
        ]),
        #Heliostat Three
        np.array([
            np.array([heliostatMirrorWidth/2+helioThreeOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([heliostatMirrorWidth/2+helioThreeOffset, heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioThreeOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioThreeOffset, heliostatMirrorHeight/2, 0.0, 1]),
        ]),
        #Heliostat Four
        np.array([
            np.array([heliostatMirrorWidth/2+helioFourOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([heliostatMirrorWidth/2+helioFourOffset, heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioFourOffset, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2+helioFourOffset, heliostatMirrorHeight/2, 0.0, 1]),
        ])
        ])

        self.helioEnclosure = np.array([
            np.array([helioEnclosureX1+helioTwoOffset, -helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX1+helioTwoOffset, helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX2+helioTwoOffset, -helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX2+helioTwoOffset, helioEnclosureY, 0.0, 1])
        ])

        self.helioNormal = np.array([0, 0, 1, 1])

    #Output: Rotates the heliostat radially (in a circle) around the power tower
    def radiallyRotate(self):
        for i in range(0, len(self.helioCoord)):
            self.helioCoord[i] = rotateCornersAroundPowerTower(self.degrees, self.helioCoord[i])

        self.helioEnclosure = rotateCornersAroundPowerTower(self.degrees, self.helioEnclosure)
        self.helioNormal = rotateAroundAbsoluteAxis(np.array([0,0,1,1]), self.degrees, self.helioNormal)

    #Output: Draws the rotated mirror set on the screen oriented in the Direction
    # that the motor rotations from TrackingV5 output.
    def rotatedMirror(self):
        #It moves the corners of the heliostat to so the normal of the
        # moved heliostat matches with the inputted normal (rotates the heliostat by
        # moving the corners)
        vertices = (
            #First rotated mirror
            (self.helioCoord[0][0][0], self.helioCoord[0][0][1], self.helioCoord[0][0][2]),
            (self.helioCoord[0][1][0], self.helioCoord[0][1][1], self.helioCoord[0][1][2]),
            (self.helioCoord[0][2][0], self.helioCoord[0][2][1], self.helioCoord[0][2][2]),
            (self.helioCoord[0][3][0], self.helioCoord[0][3][1], self.helioCoord[0][3][2]),
            #Second
            (self.helioCoord[1][0][0], self.helioCoord[1][0][1], self.helioCoord[1][0][2]),
            (self.helioCoord[1][1][0], self.helioCoord[1][1][1], self.helioCoord[1][1][2]),
            (self.helioCoord[1][2][0], self.helioCoord[1][2][1], self.helioCoord[1][2][2]),
            (self.helioCoord[1][3][0], self.helioCoord[1][3][1], self.helioCoord[1][3][2]),
            #Tself.helioCoordird
            (self.helioCoord[2][0][0], self.helioCoord[2][0][1], self.helioCoord[2][0][2]),
            (self.helioCoord[2][1][0], self.helioCoord[2][1][1], self.helioCoord[2][1][2]),
            (self.helioCoord[2][2][0], self.helioCoord[2][2][1], self.helioCoord[2][2][2]),
            (self.helioCoord[2][3][0], self.helioCoord[2][3][1], self.helioCoord[2][3][2]),
            #Fourtself.helioCoord
            (self.helioCoord[3][0][0], self.helioCoord[3][0][1], self.helioCoord[3][0][2]),
            (self.helioCoord[3][1][0], self.helioCoord[3][1][1], self.helioCoord[3][1][2]),
            (self.helioCoord[3][2][0], self.helioCoord[3][2][1], self.helioCoord[3][2][2]),
            (self.helioCoord[3][3][0], self.helioCoord[3][3][1], self.helioCoord[3][3][2]),
        )

        edges = (
            #1st helio
            (0,2),
            (2,3),
            (3,1),
            (1,0),
            #2nd
            (4,6),
            (6,7),
            (5,7),
            (4,5),
            #3rd
            (8,10),
            (10,11),
            (9,11),
            (8,9),
            #4th
            (12,14),
            (14,15),
            (13,15),
            (12,13),
        )
        glBegin(GL_QUADS)
        for edge in edges:
            glColor3fv((0.753,0.753,0.753))
            for vertex in edge:
                #Connect all the lines together by their edges
                glVertex3fv(vertices[vertex])
        glEnd()

        #This is the enclosure around the heliostat
        vertices = (
            (self.helioEnclosure[0][0], self.helioEnclosure[0][1], self.helioEnclosure[0][2]),
            (self.helioEnclosure[1][0], self.helioEnclosure[1][1], self.helioEnclosure[1][2]),
            (self.helioEnclosure[2][0], self.helioEnclosure[2][1], self.helioEnclosure[2][2]),
            (self.helioEnclosure[3][0], self.helioEnclosure[3][1], self.helioEnclosure[3][2])
        )
        edges = (
            (0,1),
            (1,3),
            (3,2),
            (2,0),
        )
        glBegin(GL_LINES)
        glColor3fv((1,1,1))
        for edge in edges:
            for vertex in edge:
                #Connect all the lines together by their edges
                glVertex3fv(vertices[vertex])
        glEnd()

    #Input: Motor rotations (x,y local rotations) to move the heliostat
    #Output: Heliostat corners after these x,y rotations and the normal of the
    # helistat before these rotations
    def tracking(self, heliostatRotationArray):

        #Restart heliostat every time it is drawn
        helio = np.array([
            np.array([heliostatMirrorWidth/2, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([heliostatMirrorWidth/2, heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2, -heliostatMirrorHeight/2, 0.0, 1]),
            np.array([-heliostatMirrorWidth/2, heliostatMirrorHeight/2, 0.0, 1]),
        ])

        self.helioEnclosure = np.array([
            np.array([helioEnclosureX1+helioTwoOffset, -helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX1+helioTwoOffset, helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX2+helioTwoOffset, -helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX2+helioTwoOffset, helioEnclosureY, 0.0, 1])
        ])
        self.helioNormal = np.array([0, 0, 1, 1])

        # Rotate the heliostat around the local X axis
        rotationMatrix = np.array(
            [[1, 0, 0, 0],
            [0, math.cos(-heliostatRotationArray[0]), -math.sin(-heliostatRotationArray[0]), 0],
            [0, math.sin(-heliostatRotationArray[0]), math.cos(-heliostatRotationArray[0]), 0],
            [0, 0, 0, 1]]
        )
        #Need to rotate the localYaxis around the X axis, then rotate the heliostat around this localYaxis
        localYaxis = np.array([0, 1, 0, 1])
        #Rotate localYaxis around green axis
        localYaxis = np.matmul(localYaxis, rotationMatrix)

        #Update the normal of the heliostat for the local X rotation, then localYaxis rotation
        self.helioNormal = np.matmul(self.helioNormal, rotationMatrix)

        #Rotate around localYaxis
        self.helioNormal = rotateAroundAbsoluteAxis(localYaxis, heliostatRotationArray[1], self.helioNormal)
        a = self.helioNormal[1]
        b = -self.helioNormal[2]
        xD = 2 * math.atan( (a - math.sqrt(a**2 + b**2)) / (b)) - math.pi / 2
        b = self.helioNormal[0]
        a = self.helioNormal[1] * math.sin(-xD) + self.helioNormal[2] * math.cos(-xD)
        yD = -2 * math.atan( (a - math.sqrt(a**2 + b**2)) / b)


        # for each corner of the heliostat, rotate the corner around the local Y, then
        # local X axis to match the output of the tracking.py
        for i in range(0, len(self.helioCoord)):
            #Rotate enclosure around x axis
            self.helioEnclosure[i] = np.matmul(self.helioEnclosure[i], rotationMatrix)

            #Rotate around x-axis (green)
            helio[i] = np.matmul(helio[i], rotationMatrix)

            #Rotate around localYaxis
            helio[i] = rotateAroundAbsoluteAxis(localYaxis, heliostatRotationArray[1], helio[i])

        self.helioCoord = np.array([
            #Heliostat Two
            np.array([
                np.array([helio[0][0]+helioTwoOffset, helio[0][1], helio[0][2], 1]),
                np.array([helio[1][0]+helioTwoOffset, helio[1][1], helio[1][2], 1]),
                np.array([helio[2][0]+helioTwoOffset, helio[2][1], helio[2][2], 1]),
                np.array([helio[3][0]+helioTwoOffset, helio[3][1], helio[3][2], 1]),
            ]),
            #Heliostat One
            np.array([
                np.array([helio[0][0]+helioOneOffset, helio[0][1], helio[0][2], 1]),
                np.array([helio[1][0]+helioOneOffset, helio[1][1], helio[1][2], 1]),
                np.array([helio[2][0]+helioOneOffset, helio[2][1], helio[2][2], 1]),
                np.array([helio[3][0]+helioOneOffset, helio[3][1], helio[3][2], 1]),
            ]),
            #Heliostat Three
            np.array([
                np.array([helio[0][0]+helioThreeOffset, helio[0][1], helio[0][2], 1]),
                np.array([helio[1][0]+helioThreeOffset, helio[1][1], helio[1][2], 1]),
                np.array([helio[2][0]+helioThreeOffset, helio[2][1], helio[2][2], 1]),
                np.array([helio[3][0]+helioThreeOffset, helio[3][1], helio[3][2], 1]),
            ]),
            #Heliostat Four
            np.array([
                np.array([helio[0][0]+helioFourOffset, helio[0][1], helio[0][2], 1]),
                np.array([helio[1][0]+helioFourOffset, helio[1][1], helio[1][2], 1]),
                np.array([helio[2][0]+helioFourOffset, helio[2][1], helio[2][2], 1]),
                np.array([helio[3][0]+helioFourOffset, helio[3][1], helio[3][2], 1]),
            ])
        ])

        #rotated helio normal before radial rotation
        copyNormal = np.array([self.helioNormal[0], self.helioNormal[1], self.helioNormal[2], 1])

        #Create the cube that is drawn based on the heliostatX/Y vectors with the enclosure at X angle
        self.radiallyRotate()
        self.rotatedMirror()

        return helio, copyNormal

#Input: point (2D array of the corners of a surface) defined by the absolute axis
#       defined by the heliostat, degrees of rotation radially around the power tower
#Output: The point rotated radially degrees around the power tower
def rotateCornersAroundPowerTower(degrees, corners):
    #Find the vector from a new coordinate system defined at the base of
    #The power tower to the point
    radians = degrees * math.pi / 180
    for i in range(0,len(corners)):
        corners[i][1] = corners[i][1] - powerTowerYoffset
        rotated = np.array([
            corners[i][0]*math.cos(radians)-corners[i][1]*math.sin(radians),
            corners[i][1]*math.cos(radians)+corners[i][0]*math.sin(radians),
            corners[i][2]*(1-math.cos(radians))+corners[i][2]*math.cos(radians),
            1
        ])
        corners[i] = rotated
        corners[i][1] = corners[i][1] + powerTowerYoffset
    return corners

#Input: radians of rotation, a point on the local axis (Green, Red, Blue)
#Output: The point (x,y,z) rotated around the power tower
def rotatePointAroundPowerTower(radians, point):
    #Find the vector from the absolute coordinate system defined at the base of
    #The power tower to the point
    temp = np.array([point[0], point[1] - powerTowerYoffset, point[2]])
    #Rotate the point around the power tower
    temp = np.array([
        temp[0]*math.cos(radians)-temp[1]*math.sin(radians),
        temp[1]*math.cos(radians)+temp[0]*math.sin(radians),
        temp[2]*(1-math.cos(radians))+temp[2]*math.cos(radians),
        1
    ])
    #Move the point back to the local coordinate system
    temp[1] = temp[1] + powerTowerYoffset
    return temp

#Input: Motor rotations for the x,y rotated heliostat
#Output: The corners of the stationary mirror after the rotations
def getStationaryMirrorCoord(localXrotation, localYrotation):

    cornerArray = [
        [stationaryMirrorWidth/2, -stationaryMirrorHeight/2, 0, 1],
        [stationaryMirrorWidth/2, stationaryMirrorHeight/2, 0, 1],
        [-stationaryMirrorWidth/2, -stationaryMirrorHeight/2, 0, 1],
        [-stationaryMirrorWidth/2, stationaryMirrorHeight/2, 0, 1]
    ]

    normal = np.array([0, 0, 1, 1])

    #Change to radians from degrees
    localXrotation = localXrotation * math.pi / 180
    localYrotation = localYrotation * math.pi / 180

    #Rotate around local X
    rotationMatrixX = np.array(
        [[1, 0, 0, 0],
        [0, math.cos(-localXrotation), -math.sin(-localXrotation), 0],
        [0, math.sin(-localXrotation), math.cos(-localXrotation), 0],
        [0, 0, 0, 1]]
    )

    rotationMatrixY = np.array(
        [[math.cos(localYrotation), 0, math.sin(localYrotation), 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [-math.sin(localYrotation), 0, math.cos(localYrotation), 0],
        [0.0, 0.0, 0.0, 1.0]]
    )

    #Rotate the mirror around the local X, then local Z axis
    for i in range(len(cornerArray)):
        cornerArray[i] = np.matmul(cornerArray[i], rotationMatrixY)
        cornerArray[i] = np.matmul(cornerArray[i], rotationMatrixX)
        #cornerArray[i] = rotateCornersAroundPowerTower(cornerArray[i], absoluteZrotation)
    normal = np.matmul(normal, rotationMatrixY)
    normal = np.matmul(normal, rotationMatrixX)

    return cornerArray, normal

#stationary mirror class holds the coordinates of each of the corners of the stationary
#mirrors.  It also has a function that will display each of the stationary mirrors
class stationaryMirror:

    #Coordinates of the stationary mirrors, Local X, Local Y
    MIRROR_ONE_ROTATION = (78.12435077, -3.37258665)
    MIRROR_TWO_ROTATION = (77.99124804, -1.377714826)
    MIRROR_THREE_ROTATION = (77.99124804, 1.377714826)
    MIRROR_FOUR_ROTATION = (78.12435077, 3.37258665)

    def __init__(self, degrees):
        self.degrees = degrees
        #Position of the stationary mirrors on the local axis
        #Temporarily have to manually type normal for the stationary mirrors - DELETE
        self.m1, self.mOneNormal = getStationaryMirrorCoord(self.MIRROR_ONE_ROTATION[0], self.MIRROR_ONE_ROTATION[1])
        self.m2, self.mTwoNormal = getStationaryMirrorCoord(self.MIRROR_TWO_ROTATION[0], self.MIRROR_TWO_ROTATION[1])
        self.m3, self.mThreeNormal = getStationaryMirrorCoord(self.MIRROR_THREE_ROTATION[0], self.MIRROR_THREE_ROTATION[1])
        self.m4, self.mFourNormal = getStationaryMirrorCoord(self.MIRROR_FOUR_ROTATION[0], self.MIRROR_FOUR_ROTATION[1])
        #Add distance from helio to stationary to y
        #And add distance in the x axis - offset between mirrors
        for row in range(len(self.m1)):
            self.m1[row][1] += distanceStationaryMirrorToRotatedMirror
            self.m1[row][0] += helioOneOffset
            self.m2[row][1] += distanceStationaryMirrorToRotatedMirror
            self.m2[row][0] += helioTwoOffset
            self.m3[row][1] += distanceStationaryMirrorToRotatedMirror
            self.m3[row][0] += helioThreeOffset
            self.m4[row][1] += distanceStationaryMirrorToRotatedMirror
            self.m4[row][0] += helioFourOffset

        #Initialize the enclosure around the mirrors
        #Enclosure around the stationary mirrors
        self.enclosure = np.array([
            np.array([helioEnclosureX1+helioTwoOffset, -helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX1+helioTwoOffset, helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX2+helioTwoOffset, -helioEnclosureY, 0.0, 1]),
            np.array([helioEnclosureX2+helioTwoOffset, helioEnclosureY, 0.0, 1])
        ])

        angleToMoveRadians = self.MIRROR_ONE_ROTATION[0] * math.pi / 180
        rotationMatrix = np.array(
            [[1, 0, 0, 0],
            [0, math.cos(-angleToMoveRadians), -math.sin(-angleToMoveRadians), 0],
            [0, math.sin(-angleToMoveRadians), math.cos(-angleToMoveRadians), 0],
            [0, 0, 0, 1]]
        )
        for i in range(0, len(self.enclosure)):
            #Rotate enclosure around x axis
            self.enclosure[i] = np.matmul(self.enclosure[i], rotationMatrix)
            self.enclosure[i][1] = self.enclosure[i][1] + distanceStationaryMirrorToRotatedMirror

    #Output: Draws the stationary mirrors rotated radially around the power tower
    def radiallyRotate(self):
        self.m1 = rotateCornersAroundPowerTower(self.degrees, self.m1)
        self.m2 = rotateCornersAroundPowerTower(self.degrees, self.m2)
        self.m3 = rotateCornersAroundPowerTower(self.degrees, self.m3)
        self.m4 = rotateCornersAroundPowerTower(self.degrees, self.m4)

        self.enclosure = rotateCornersAroundPowerTower(self.degrees, self.enclosure)

    def getMirrorOne(self):
        return self.m1
    def getMirrorTwo(self):
        return self.m2
    def getMirrorThree(self):
        return self.m3
    def getMirrorFour(self):
        return self.m4
    def getNormalOne(self):
        return self.mOneNormal
    def getNormalTwo(self):
        return self.mTwoNormal
    def getNormalThree(self):
        return self.mThreeNormal
    def getNormalFour(self):
        return self.mFourNormal

    #Output: Draws the stationary mirror set
    def stationaryMirror(self):
        vertices = (
            #First stationary mirror
            (self.m2[0][0], self.m2[0][1], self.m2[0][2]),
            (self.m2[1][0], self.m2[1][1], self.m2[1][2]),
            (self.m2[2][0], self.m2[2][1], self.m2[2][2]),
            (self.m2[3][0], self.m2[3][1], self.m2[3][2]),
            #Second
            (self.m1[0][0], self.m1[0][1], self.m1[0][2]),
            (self.m1[1][0], self.m1[1][1], self.m1[1][2]),
            (self.m1[2][0], self.m1[2][1], self.m1[2][2]),
            (self.m1[3][0], self.m1[3][1], self.m1[3][2]),
            #Third
            (self.m3[0][0], self.m3[0][1], self.m3[0][2]),
            (self.m3[1][0], self.m3[1][1], self.m3[1][2]),
            (self.m3[2][0], self.m3[2][1], self.m3[2][2]),
            (self.m3[3][0], self.m3[3][1], self.m3[3][2]),
            #Fourth
            (self.m4[0][0], self.m4[0][1], self.m4[0][2]),
            (self.m4[1][0], self.m4[1][1], self.m4[1][2]),
            (self.m4[2][0], self.m4[2][1], self.m4[2][2]),
            (self.m4[3][0], self.m4[3][1], self.m4[3][2]),
        )
        edges = (
            #1st mirror
            (0,2),
            (2,3),
            (1,3),
            (0,1),
            #2nd
            (4,6),
            (6,7),
            (5,7),
            (4,5),
            #3rd
            (8,10),
            (10,11),
            (9,11),
            (8,9),
            #4th
            (12,14),
            (14,15),
            (13,15),
            (12,13),
        )
        glBegin(GL_QUADS)
        for edge in edges:
            glColor3fv((0.753,0.753,0.753))
            for vertex in edge:
                #Connect all the lines together by their edges
                glVertex3fv(vertices[vertex])
        glEnd()

        #This is the self.enclosure around the heliostat
        vertices = (
            (self.enclosure[0][0], self.enclosure[0][1], self.enclosure[0][2]),
            (self.enclosure[1][0], self.enclosure[1][1], self.enclosure[1][2]),
            (self.enclosure[2][0], self.enclosure[2][1], self.enclosure[2][2]),
            (self.enclosure[3][0], self.enclosure[3][1], self.enclosure[3][2])
        )
        edges = (
            (0,2),
            (0,1),
            (1,3),
            (2,3),
        )
        glBegin(GL_LINES)
        glColor3fv((1,1,1))
        for edge in edges:
            for vertex in edge:
                #Connect all the lines together by their edges
                glVertex3fv(vertices[vertex])
        glEnd()

#Rays are simulated light beams from the sun that bounce off the mirrors.
class Ray:

    def __init__(self, direction, degree):
        self.reflected = False
        self.direction = direction
        self.degree = degree
        self.nextReflected = None
        self.color = (0.8, 0.8, 0)
        self.length = 1000.0
        self.startingLocation = [0.0, 0.0, 0.0]

    def isReflected(self):
        return self.reflected

    def setReflected(self, bool):
        self.reflected = bool

    def getDirection(self):
        return self.direction

    def setDirection(self, vertex):
        self.direction = vertex

    def getNextReflected(self):
        return self.nextReflected

    def setNextReflected(self, Ray):
        self.nextReflected = Ray

    def getLength(self):
        return self.length

    def setLength(self, length):
        self.length = length

    def getStartingLocation(self):
        return self.startingLocation

    def setStartingLocation(self, startingLocation):
        self.startingLocation = startingLocation

    #Rotate the light rays around the power tower specified by the global
    # degrees variable
    def radiallyRotate(self):
        radians = self.degree * math.pi / 180
        self.setStartingLocation(rotatePointAroundPowerTower(radians, self.startingLocation))
        self.nextReflected.startingLocation = rotatePointAroundPowerTower(radians, self.nextReflected.startingLocation)
        self.nextReflected.direction = rotateAroundAbsoluteAxis(np.array([0,0,1,1]), radians, self.nextReflected.direction)

        if self.nextReflected.nextReflected is not None:
            self.nextReflected.nextReflected.startingLocation = rotatePointAroundPowerTower(radians, self.nextReflected.nextReflected.startingLocation)
            self.nextReflected.nextReflected.direction = rotateAroundAbsoluteAxis(np.array([0,0,1,1]), radians, self.nextReflected.nextReflected.direction)
    #Draws all of the vectors that are in this linked list
    def draw(self):

        line = (
            (self.startingLocation[0], self.startingLocation[1], self.startingLocation[2]),
            (self.direction[0] * self.length + self.startingLocation[0], self.direction[1] * self.length + self.startingLocation[1], self.direction[2] * self.length + self.startingLocation[2])
        )
        lineEdge = (
            (0,1),
        )
        glBegin(GL_LINES)
        glColor3fv(self.color)
        for edge in lineEdge:
            for vertex in edge:
                glVertex3fv(line[vertex])
        glEnd()

#Input: Vector that points from the local coordinate system
# to the sun
#Output: Displays a 'sun' at the inputted vector location
def sun(sunRotation):

    #The only way to create a sphere is to have this shit
    sphere = gluNewQuadric()

    distanceAway = 1000

    glPushMatrix()

    #Move to where the sun will be located
    glTranslatef(sunRotation[0]*distanceAway, sunRotation[1]*distanceAway, sunRotation[2]*distanceAway) #Move to the place
    glColor4f(0.9, 0.2, 0.2, 1) #Put color
    gluSphere(sphere, 30, 32, 16) #Draw sphere

    glPopMatrix()

#Draws the power tower
def powerTower():

    global powerTowerYoffset
    global powerTowerZoffset
    global powerTowerDimension

    #Light catcher w/ dimensions powerTowerDimension ^ 3
    vertices = (
        (powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, powerTowerZoffset-powerTowerDimension/2),
        (powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset-powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset-powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, powerTowerZoffset-powerTowerDimension/2),
        (powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, powerTowerZoffset+powerTowerDimension/2),
        (powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset+powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, powerTowerZoffset+powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset+powerTowerDimension/2)
    )
    surfaces = (
        (0,1,2,3),
        (0,1,5,4),
        (0,3,6,4),
        (4,5,7,6),
        (1,2,7,5),
        (3,2,7,6)
    )
    glBegin(GL_QUADS)
    for surface in surfaces:
         glColor3fv((0, 0.5, 0.5))
         for vertex in surface:
             glVertex3fv(vertices[vertex])
    glEnd()

    #The power tower holder
    vertices = (
        (powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, powerTowerZoffset-powerTowerDimension/2),
        (powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset-powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset-powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, powerTowerZoffset-powerTowerDimension/2),
        (powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, -powerTowerDimension/2),
        (powerTowerDimension/2, powerTowerYoffset, -powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset-powerTowerDimension, -powerTowerDimension/2),
        (-powerTowerDimension/2, powerTowerYoffset, -powerTowerDimension/2)
    )
    surfaces = (
        (0,1,2,3),
        (0,1,5,4),
        (0,3,6,4),
        (4,5,7,6),
        (1,2,7,5),
        (3,2,7,6)
    )
    glBegin(GL_QUADS)
    for surface in surfaces:
         glColor3fv((0.754,0.754,0.754))
         for vertex in surface:
             glVertex3fv(vertices[vertex])
    glEnd()

#Input:  Axis:
#Returns 'rotatedVector' is rotated 'radians' radians around 'axis' vector
def rotateAroundAbsoluteAxis(axis, radians, rotatedVector):
    #Matrix to rotate the rotatedVector around the axis
    rotationMatrix = np.array(
        [[((1-math.cos(radians))*axis[0]**2+math.cos(radians)), ((1-math.cos(radians))*axis[0]*axis[1]-math.sin(radians)*axis[2]), ((1-math.cos(radians))*axis[0]*axis[2]+math.sin(radians)*axis[1]), 0],
        [((1-math.cos(radians))*axis[0]*axis[1]+math.sin(radians)*axis[2]), ((1-math.cos(radians))*axis[1]**2+math.cos(radians)), ((1-math.cos(radians))*axis[1]*axis[2]-math.sin(radians)*axis[0]), 0],
        [((1-math.cos(radians))*axis[0]*axis[2]-math.sin(radians)*axis[1]), ((1-math.cos(radians))*axis[1]*axis[2]+math.sin(radians)*axis[0]), ((1-math.cos(radians))*axis[2]**2)+math.cos(radians), 0],
        [0, 0, 0, 1]]
    )
    return np.matmul(rotationMatrix, rotatedVector)

#offset in hours of the current time in Madison, WI
#In case it is late you can push back the time and have live sunlight
offset = 0 #You need to add a constant to get to current time. California constant = 9.85

#What degrees should the heliostat system be rotated radially around the power tower
degrees = 90 * math.pi / 180

# Calls the trackingV5 algorithm's main and returns the array of absolute visual x
# rotation degrees, absolute visual z rotation degrees
#
# Returns the angles of rotation in radians
def getTrackingCoordinates():
    LST=90 #Standard Meridian
    Latitude=33.590395
    Longitude=-114.53423
    fp='./MirrorInformationV4.csv'

    global offset
    global degrees
    tOffset = offset
    tDegrees = degrees
    sunRotation, solarTime = trackingV5.main(LST,Longitude,Latitude,fp, tOffset)
    sunRotation = rotateAroundAbsoluteAxis(np.array([0,0,1,1]), tDegrees, sunRotation)
    normal = getNormalVector(sunRotation)
    xRotation, yRotation = getRotationAngles(normal)
    #sunRotation = rotateAroundAbsoluteAxis(np.array([0,0,1,1]), -tDegrees, sunRotation)

    #print("Tracking array: "+str(trackingArray))
    #print("Sun rotation: "+str(sunRotation))
    #This is just to make out codes work together, trackingVX defines sun as
    #sun to origin, this defines as origin to sun
    sunRotation[0:3] = sunRotation[0:3] * -1

    offset = offset + 0.01
    return [xRotation, yRotation], sunRotation, solarTime
    #return [float(trackingArray[0][1]) * math.pi / 180, float(trackingArray[0][2]) * math.pi / 180], sunRotation

#Input: Incoming vector hitting a plane defined by the normal vector
#Output: Reflected vector of the incoming vector off a plane defined by the normal vector
def getReflectedVector(incomingVector, normal):
    #Find the reflected angle of the first beams - the angle reflected is the same for all beams
    #Law of reflection: R = 2(N dot L)N - L
    dotProduct = np.dot(normal[0:3], (incomingVector[0:3]))
    mult = 2 * dotProduct * normal[0:3]
    #Need to flip direction because this output is facing the wrong direction.
    reflectedVector = np.subtract(mult, (incomingVector[0:3])) * -1
    reflectedVector = np.array([reflectedVector[0], reflectedVector[1], reflectedVector[2], 1])
    return reflectedVector

#Input: Incoming vector, with outgoing vector pointing in the y direction
#Output: Normal from the plane that would reflect the incoming vector to
# the outgoing vector (in y-axis)
def getNormalVector(incomingVector):
    incoming = -1 * np.array([incomingVector[0], incomingVector[1], incomingVector[2]])
    outgoingVector = np.array([0,1,0])
    mult = np.add(outgoingVector, incoming)

    dotSqr = np.dot(mult/2, incoming)
    dot = math.sqrt(dotSqr)

    normal = mult / (2 * dot)
    normal[2] = abs(normal[2])
    return normal

#From a normal vector, find the x and y local rotations necessary to point in the
# normal direction
def getRotationAngles(normal):
    a = normal[1]
    b = -normal[2]
    xD = 2 * math.atan( (a - math.sqrt(a**2 + b**2)) / (b)) - math.pi / 2

    b = normal[0]
    a = normal[1] * math.sin(-xD) + normal[2] * math.cos(-xD)
    yD = -2 * math.atan( (a - math.sqrt(a**2 + b**2)) / b)

    return xD, yD

#This is used to find where the light rays should start on the rotated heliostat
#Input: fraction of the heliostat (%) where the ray should hit the rotated mirror,
# and the normal of the helistat mirror plane, with a point on the plane in order
# to fully define the plane
#Output: Coordinates (x,y,z) where the ray should start based on the local coordinate
# system
def rayStartingLocation(fractionHeight, fractionWidth, normal, helioCorners):
    fractionHeightScaled = (float(fractionHeight) - 0.5) * 2
    fractionWidthScaled = (float(fractionWidth) - 0.5) * 2
    x = fractionWidthScaled * helioCorners[0][0]
    y = fractionHeightScaled * helioCorners[0][1]
    #Need to shift down due to local Y rotation
    tmpX = x
    tmpY = y
    ratio = abs(helioCorners[0][1] / helioCorners[0][0])
    if(helioCorners[3][1] > helioCorners[1][1]):
        #If along the right side of the diagonal of the heliostat
        if(x * ratio + y > 0):
            #Need offset for the y rotation
            dy = (helioCorners[3][1] - helioCorners[1][1]) * fractionWidth * (1 - fractionHeight)
            y -= dy
        #If along the left side of the diagonal of the heliostat
        elif(tmpX * ratio + tmpY < 0):
            #Need offset for the y rotation
            dy = (helioCorners[3][1] - helioCorners[1][1]) * (1 - fractionWidth) * fractionHeight
            y += dy
    else:
        #If along the right side of the diagonal of the heliostat
        if(x * ratio + y > 0):
            #Need offset for the y rotation
            dy = (helioCorners[1][1] - helioCorners[3][1]) * fractionWidth * (1 - fractionHeight)
            y += dy
        #If along the left side of the diagonal of the heliostat
        elif(tmpX * ratio + tmpY < 0):
            #Need offset for the y rotation
            dy = (helioCorners[1][1] - helioCorners[3][1]) * (1 - fractionWidth) * fractionHeight
            y -= dy
    z = (-normal[0] * x - normal[1] * y) / normal[2]
    return [x, y, z]

#Input: a light ray class that may or may not hit the mirror of the stationary mirror,
# normal of the stationary mirrors, and a point on the stationary mirror
# plane.
#Output: boolean of whether or not the incoming ray hits the stationary mirror
#If the incoming ray hits the stationary mirror, it also returns the length away from the
#rotated mirror
#
#Defines a plane with the stationary mirror normal and one of the corner coordinates.
#Finds where the ray direction intersects with the plane, then finds the disance
#in the x and y direction, comparing to see if it is in range of the mirror coordinates
def stationaryMirrorRefected(ray, smCoords):
    #Part A:  Is the ray within Right edge of the mirror?
    m = (smCoords[1][0]-smCoords[0][0])/(smCoords[1][2]-smCoords[0][2])
    b = smCoords[0][0] - m * smCoords[0][2]
    if(m * ray.startingLocation[2] + b < ray.startingLocation[0]):
        return False

    #Part B: Is the ray above the bottom edge of the mirror?
    m = (smCoords[0][2]-smCoords[2][2])/(smCoords[0][0]-smCoords[2][0])
    b = smCoords[0][2] - m * smCoords[0][0]
    if(m * ray.startingLocation[0] + b > ray.startingLocation[2]):
        return False

    #Part C: Is the ray to the within the left edge of the mirror?
    m = (smCoords[3][0]-smCoords[2][0])/(smCoords[3][2]-smCoords[2][2])
    b = smCoords[2][0] - m * smCoords[2][2]
    if(m * ray.startingLocation[2] + b > ray.startingLocation[0]):
        return False

    #Part D: Is the ray below the tom edge?
    m = (smCoords[1][2]-smCoords[3][2])/(smCoords[1][0]-smCoords[3][0])
    b = smCoords[1][2] - m * smCoords[1][0]
    if(m * ray.startingLocation[0] + b < ray.startingLocation[2]):
        return False
    return True

#Input: a light ray class, and a points on a plane, and the normal of the plane
#Output: Find distance between the plane and the starting location of the ray
def lengthToPlane(ray, plane, normal):
    #Define a plane (ax+by+cz = d)
    #a,b,c are normal coordinates
    #Use the direction vector (x = 1 + 2*t) for x,y,z
    #Solve for t
    #Plug in t for x,y,z

    #Starting Location
    sL = ray.getStartingLocation()
    #Direction
    dir = ray.getDirection()
    d = normal[0]*plane[0][0] + normal[1]*plane[0][1] + normal[2]*plane[0][2]
    #Plug in substitution for x,y,z (eg x=1+2*t)
    tSubstitution = normal[0]*sL[0] + normal[1]*sL[1] + normal[2]*sL[2]
    #Use vector direction to find constant
    vectorConstant = normal[0]*dir[0] + normal[1]*dir[1] + normal[2]*dir[2]
    t = (d-tSubstitution)/vectorConstant
    #These are the coordinates of the location where the light lands on the stationary mirror
    x = sL[0] + dir[0]*t
    y = sL[1] + dir[1]*t
    z = sL[2] + dir[2]*t

    #Find the distance between the starting location and the end location
    length = math.sqrt((sL[0] - x)**2 + (sL[1] - y)**2 + (sL[2] - z)**2)
    return float(length), [x,y,z]

#Draws the inputted ray
def drawRay(currentRay):
    if currentRay is None:
        return
    else:
        currentRay.draw()
        drawRay(currentRay.getNextReflected())

# Game class contains the actual pygame environment, all the local/absolute coordinate systems
# for the heliostat, and all of the matrix transformations.
#
# MouseMove:  Allows for the mouse to click and drag the screen around the absolute visual coordinate system
# Main:  Calls all of the objects (rectangle, coordinate system, etc) that are drawn, and updates the drawing
#        indefinetly.
class Game(object):

    def __init__(self):

        global degrees
        tDegrees = degrees

        # Unit vectors that are used for visual turning of the screen
        # They go from the absolute visual vector to the local vector
        self.xAbsoluteToLocal = np.array([1,0,0,1])
        self.yAbsoluteToLocal = np.array([0,1,0,1])

        self.degrees = 0
        #Create one heliostat at a rotation of zero degrees around the power tower
        self.rotatedMirror = rotatedMirror(self.degrees)

        #Used in the ray tracing, needs coordinates of each of the mirrors
        self.stationaryMirror = stationaryMirror(self.degrees)

        #Vector that represents the angle of the sun compared to the local axis
        self.sunRotation = None

        self.solarTime = None

        #Must be int
        self.numberOfRays = 5

        #This array holds all the of the heads of the linked list of Rays that start at the Heliostat
        #Is initialized when tracking (p) is turned on
        self.rayArray2D_1 = [[None for i in range(self.numberOfRays)] for j in range(self.numberOfRays)]
        self.rayArray2D_2 = [[None for i in range(self.numberOfRays)] for j in range(self.numberOfRays)]
        self.rayArray2D_3 = [[None for i in range(self.numberOfRays)] for j in range(self.numberOfRays)]
        self.rayArray2D_4 = [[None for i in range(self.numberOfRays)] for j in range(self.numberOfRays)]

        # Used to move the visual screen around with the mouse and find how much to move the screen by
        self.lastMousePositionX = 0
        self.lastMousePositionY = 0

    #Input: a matrix that rotates a vector around the x,y or z axis
    # This is used to rotate the screen when the user clicks and drags with the mouse
    def updateVectorsLocalRotation(self, rotationMatrix):
        #When the local axis is moved, the visual absolute vertex needs to be updated.
        self.yAbsoluteToLocal = np.matmul(self.yAbsoluteToLocal, rotationMatrix)
        self.xAbsoluteToLocal = np.matmul(self.xAbsoluteToLocal, rotationMatrix)

    # Used to rotate the enviroment around the local (X - Green, Y - Red, Z - Blue)
    # Axis's
    # Use q,w to rotate around the positive, negative X axis
    # Use a,s to rotate around the positive, negative Y axis
    # Use z,x to rotate around the positive, negative Z axis
    # Use e to rotate -90 around x (Use first to orient Z in the up direction (toward sky))
    def localAxisRotations(self, event):

        #In degrees
        angleToMove = 2
        angleToMoveRadians = angleToMove * math.pi / 180

        #Rotate locally around the positive y-axis
        if(event.type == pygame.KEYDOWN and event.key == pygame.K_a):

            rotationMatrix = np.array(
                [[math.cos(angleToMoveRadians), 0, math.sin(angleToMoveRadians), 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [-math.sin(angleToMoveRadians), 0, math.cos(angleToMoveRadians), 0],
                [0.0, 0.0, 0.0, 1.0]]
            )

            #Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            #Rotate locally around the positive Y
            glRotatef(angleToMove,0,1,0)

        #Rotate around the negative local y axis
        if(event.type == pygame.KEYDOWN and event.key == pygame.K_s):

            rotationMatrix = np.array(
                [[math.cos(-angleToMoveRadians), 0, math.sin(-angleToMoveRadians), 0],
                [0, 1, 0, 0],
                [-math.sin(-angleToMoveRadians), 0, math.cos(-angleToMoveRadians), 0],
                [0, 0, 0, 1]]
            )

            #Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            #Rotate around the negative local y axis
            glRotatef(-angleToMove,0,1,0)

        #Rotate around the positive local x axis
        if(event.type == pygame.KEYDOWN and event.key == pygame.K_q):

            rotationMatrix = np.array(
                [[1, 0, 0, 0],
                [0, math.cos(angleToMoveRadians), -math.sin(angleToMoveRadians), 0],
                [0, math.sin(angleToMoveRadians), math.cos(angleToMoveRadians), 0],
                [0, 0, 0, 1]]
            )

            #Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            #Rotate around the positive local x axis
            glRotatef(angleToMove,1,0,0)

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_w):

            rotationMatrix = np.array(
                [[1, 0, 0, 0],
                [0, math.cos(-angleToMoveRadians), -math.sin(-angleToMoveRadians), 0],
                [0, math.sin(-angleToMoveRadians), math.cos(-angleToMoveRadians), 0],
                [0, 0, 0, 1]]
            )

            #Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            #Rotate around the negative local x axis
            glRotatef(-angleToMove,1,0,0)

        # If key == e, rotate -90 around X, -90 around Z
        if(event.type == pygame.KEYDOWN and event.key == pygame.K_e):

            angleRotation = -90 * math.pi / 180
            rotationMatrix = np.array(
                [[1, 0, 0, 0],
                [0, math.cos(angleRotation), -math.sin(angleRotation), 0],
                [0, math.sin(angleRotation), math.cos(angleRotation), 0],
                [0, 0, 0, 1]]
            )

            # Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            # Rotate around -Z axis 90 degrees
            rotationMatrix = np.array(
                [[math.cos(angleRotation), -math.sin(angleRotation), 0, 0],
                [math.sin(angleRotation), math.cos(angleRotation), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]
            )

            self.updateVectorsLocalRotation(rotationMatrix)

            # Rotate around the negative local x axis
            glRotatef(-90,1,0,0)
            # Rotate around the negatie local Z axis
            glRotatef(-90,0,0,1)

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_z):

            rotationMatrix = np.array(
                [[math.cos(angleToMoveRadians), -math.sin(angleToMoveRadians), 0, 0],
                [math.sin(angleToMoveRadians), math.cos(angleToMoveRadians), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]
            )

            #Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            #Rotate around the positive local z axis
            glRotatef(angleToMove,0,0,1)

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_x):

            rotationMatrix = np.array(
                [[math.cos(-angleToMoveRadians), -math.sin(-angleToMoveRadians), 0, 0],
                [math.sin(-angleToMoveRadians), math.cos(-angleToMoveRadians), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]]
            )

            #Change all self. vectors for rotation
            self.updateVectorsLocalRotation(rotationMatrix)

            #Rotate around the negative local Z axis
            glRotatef(-angleToMove,0,0,1)

    # Used to rotate the environment using the arrow keys around the absolute visual
    # Coordinate system
    def absoluteRotate(self, event):

        #In degrees
        angleToMove = 1
        angleToMoveRadians = angleToMove * math.pi / 180

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_DOWN):

            #Rotate around the negative X
            glRotatef(angleToMove, self.xAbsoluteToLocal[0], self.xAbsoluteToLocal[1], self.xAbsoluteToLocal[2])

            #Matrix to rotate the absolute visual y axis around the negative absolute visual x axis
            self.yAbsoluteToLocal = rotateAroundAbsoluteAxis(self.xAbsoluteToLocal, -angleToMoveRadians, self.yAbsoluteToLocal)

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_UP):

            #Rotate around the positive X
            glRotatef(-angleToMove, self.xAbsoluteToLocal[0], self.xAbsoluteToLocal[1], self.xAbsoluteToLocal[2])

            #Matrix to rotate the absolute visual y axis around the positive absolute visual X axis
            self.yAbsoluteToLocal = rotateAroundAbsoluteAxis(self.xAbsoluteToLocal, angleToMoveRadians, self.yAbsoluteToLocal)

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_RIGHT):

            #Rotate around the positive Y axis
            glRotatef(-angleToMove, self.yAbsoluteToLocal[0], self.yAbsoluteToLocal[1], self.yAbsoluteToLocal[2])

            #Matrix to rotate the absolute visual x axis around the positive absolute visual y axis
            self.xAbsoluteToLocal = rotateAroundAbsoluteAxis(self.yAbsoluteToLocal, angleToMoveRadians, self.xAbsoluteToLocal)

        if(event.type == pygame.KEYDOWN and event.key == pygame.K_LEFT):

            #Rotate around the negative Y axis
            glRotatef(angleToMove, self.yAbsoluteToLocal[0], self.yAbsoluteToLocal[1], self.yAbsoluteToLocal[2])

            #Matrix to rotate the absolute visual x axis around the negative absolute visual y axis
            self.xAbsoluteToLocal = rotateAroundAbsoluteAxis(self.yAbsoluteToLocal, -angleToMoveRadians, self.xAbsoluteToLocal)

    # Can use mouse + click + drag to rotate the entire environment
    def mouseMove(self, event):

        #In degrees
        angleToMove = 5
        angleToMoveRadians = angleToMove * math.pi / 180

        if(event.type == pygame.MOUSEMOTION):
            x, y = pygame.mouse.get_pos()
            deltaX = self.lastMousePositionX - x
            deltaY = self.lastMousePositionY - y

            #If mouse if clicked
            if((pygame.mouse.get_pressed())[0]):
                #Rotate visual coordinates up
                if(deltaX > 0):
                    #Rotate around the positive Y axis
                    glRotatef(-angleToMove, self.yAbsoluteToLocal[0], self.yAbsoluteToLocal[1], self.yAbsoluteToLocal[2])

                    #Matrix to rotate the absolute visual x axis around the positive absolute visual y axis
                    self.xAbsoluteToLocal = rotateAroundAbsoluteAxis(self.yAbsoluteToLocal, angleToMoveRadians, self.xAbsoluteToLocal)

                if(deltaX < 0):
                    #Rotate around the negative Y axis
                    glRotatef(angleToMove, self.yAbsoluteToLocal[0], self.yAbsoluteToLocal[1], self.yAbsoluteToLocal[2])

                    #Matrix to rotate the absolute visual x axis around the negative absolute visual y axis
                    self.xAbsoluteToLocal = rotateAroundAbsoluteAxis(self.yAbsoluteToLocal, -angleToMoveRadians, self.xAbsoluteToLocal)

                if(deltaY > 0):
                    #Rotate around the positive X
                    glRotatef(-angleToMove, self.xAbsoluteToLocal[0], self.xAbsoluteToLocal[1], self.xAbsoluteToLocal[2])

                    #Matrix to rotate the absolute visual y axis around the positive absolute visual X axis
                    self.yAbsoluteToLocal = rotateAroundAbsoluteAxis(self.xAbsoluteToLocal, angleToMoveRadians, self.yAbsoluteToLocal)

                if(deltaY < 0):
                    #Rotate around the negative X
                    glRotatef(angleToMove, self.xAbsoluteToLocal[0], self.xAbsoluteToLocal[1], self.xAbsoluteToLocal[2])

                    #Matrix to rotate the absolute visual y axis around the negative absolute visual x axis
                    self.yAbsoluteToLocal = rotateAroundAbsoluteAxis(self.xAbsoluteToLocal, -angleToMoveRadians, self.yAbsoluteToLocal)

            self.lastMousePositionX = x
            self.lastMousePositionY = y

    # Use the -,= keys to zoom in, out
    def zoomInOut(self, event):
        #Scale in and out of the environment
        if(event.type == pygame.KEYDOWN and event.key == pygame.K_EQUALS):
            glScaled(1.1, 1.1, 1.1)
        if(event.type == pygame.KEYDOWN and event.key == pygame.K_MINUS):
            glScaled(0.9, 0.9, 0.9)

    #Input: rotated heliostat mirror corner locations, normal vector of rotated mirror normal
    #Output: Draws all light lines that bounce from the sun to rotated mirror to
    # stationary mirror (if the light ray hits the mirror)
    def rayTracing(self, rotatedHelioCoord, rotatedNormal):

        #Unit vector of reflection angle off rotated helio, same for all Rays
        helioReflection = getReflectedVector(self.sunRotation[0:3] * -1, self.rotatedMirror.helioNormal)
        stationaryReflection1 = getReflectedVector(helioReflection, self.stationaryMirror.getNormalOne())
        stationaryReflection2 = getReflectedVector(helioReflection, self.stationaryMirror.getNormalTwo())
        stationaryReflection3 = getReflectedVector(helioReflection, self.stationaryMirror.getNormalThree())
        stationaryReflection4 = getReflectedVector(helioReflection, self.stationaryMirror.getNormalFour())

        powerTowerSurface = [[powerTowerDimension/2, powerTowerYoffset, powerTowerZoffset+powerTowerDimension/2]]

        #Initialize the head of the linked list of Ray with the sunRotation vector
        for row in range(len(self.rayArray2D_2)):
            for col in range(len(self.rayArray2D_2[row])):
                #The starting location is also the same for all the heliostats, only requiring an offset
                startingLocation = rayStartingLocation(float(row)/(len(self.rayArray2D_2)-1), float(col)/(len(self.rayArray2D_2[row])-1), rotatedNormal, rotatedHelioCoord)
                #Create ray tracing for heliostat 2
                self.rayArray2D_2[row][col] = Ray(self.sunRotation, self.degrees)
                #Offset the starting location for each mirror
                startingTwo = np.array([startingLocation[0] + helioTwoOffset, startingLocation[1], startingLocation[2], 1])
                self.rayArray2D_2[row][col].setStartingLocation(startingTwo)
                #Set next ray in the linked list to the reflection of the sun's ray off the rotated heliostat
                self.rayArray2D_2[row][col].setNextReflected(Ray(helioReflection, self.degrees))

                #Set the next ray starting location (it is the same as the sun ray)
                self.rayArray2D_2[row][col].getNextReflected().setStartingLocation(startingTwo)
                #Set next ray distance of heliostat to stationary mirror
                length, stationStartingCoord = lengthToPlane(self.rayArray2D_2[row][col].getNextReflected(), self.stationaryMirror.getMirrorTwo(), self.stationaryMirror.getNormalTwo())
                #Set the next ray length
                self.rayArray2D_2[row][col].getNextReflected().setLength(length)
                #Create ray between stationary mirror and power tower
                # IF the ray hits the stationary mirrors
                if(stationaryMirrorRefected(self.rayArray2D_2[row][col].getNextReflected(), self.stationaryMirror.m2)):
                    rayBetweenStationPower = Ray(stationaryReflection2, self.degrees)
                    self.rayArray2D_2[row][col].getNextReflected().setNextReflected(rayBetweenStationPower)
                    self.rayArray2D_2[row][col].getNextReflected().reflected = True
                    rayBetweenStationPower.setStartingLocation(stationStartingCoord)
                    rayBetweenStationPower.setLength(lengthToPlane(rayBetweenStationPower, powerTowerSurface, [0,1,0,1])[0])
                else:
                    self.rayArray2D_2[row][col].getNextReflected().setLength(length*100)

                self.rayArray2D_2[row][col].radiallyRotate()
                drawRay(self.rayArray2D_2[row][col])

                #Create mirror three and set parameters for light from sun to heliostat
                self.rayArray2D_3[row][col] = Ray(self.sunRotation, self.degrees)
                startingThree = np.array([startingLocation[0] + helioThreeOffset, startingLocation[1], startingLocation[2], 1])
                #Set the starting location on the rotated heliostat for every point
                self.rayArray2D_3[row][col].setStartingLocation(startingThree)
                #Set next ray in the linked list to the reflection of the sun's ray off the rotated heliostat
                self.rayArray2D_3[row][col].setNextReflected(Ray(helioReflection, self.degrees))

                #Set the next ray starting location (it is the same as the sun ray)
                self.rayArray2D_3[row][col].getNextReflected().setStartingLocation(startingThree)
                #Set next ray distance of heliostat to stationary mirror
                length, stationStartingCoord = lengthToPlane(self.rayArray2D_3[row][col].getNextReflected(), self.stationaryMirror.getMirrorThree(), self.stationaryMirror.getNormalThree())
                #Set the next ray length
                self.rayArray2D_3[row][col].getNextReflected().setLength(length)

                #Check if the ray actually hits the stationary mirror
                if(stationaryMirrorRefected(self.rayArray2D_3[row][col].getNextReflected(), self.stationaryMirror.m3)):
                    rayBetweenStationPower = Ray(stationaryReflection3, self.degrees)
                    self.rayArray2D_3[row][col].getNextReflected().setNextReflected(rayBetweenStationPower)
                    self.rayArray2D_3[row][col].getNextReflected().reflected = True
                    rayBetweenStationPower.setStartingLocation(stationStartingCoord)
                    rayBetweenStationPower.setLength(lengthToPlane(rayBetweenStationPower, powerTowerSurface, [0,1,0,1])[0])
                else:
                    self.rayArray2D_3[row][col].getNextReflected().setLength(length*100)
                self.rayArray2D_3[row][col].radiallyRotate()
                drawRay(self.rayArray2D_3[row][col])


                self.rayArray2D_1[row][col] = Ray(self.sunRotation, self.degrees)
                startingOne = np.array([startingLocation[0] + helioOneOffset, startingLocation[1], startingLocation[2], 1])
                #Set the starting location on the rotated heliostat for every point
                self.rayArray2D_1[row][col].setStartingLocation(startingOne)
                #Set next ray in the linked list to the reflection of the sun's ray off the rotated heliostat
                self.rayArray2D_1[row][col].setNextReflected(Ray(helioReflection, self.degrees))

                #Set the next ray starting location (it is the same as the sun ray)
                self.rayArray2D_1[row][col].getNextReflected().setStartingLocation(startingOne)
                #Set next ray distance of heliostat to stationary mirror
                length, stationStartingCoord = lengthToPlane(self.rayArray2D_1[row][col].getNextReflected(), self.stationaryMirror.getMirrorOne(), self.stationaryMirror.getNormalOne())
                #Set the next ray length
                self.rayArray2D_1[row][col].getNextReflected().setLength(length)

                #Check if the ray actually hits the stationary mirror
                if(stationaryMirrorRefected(self.rayArray2D_1[row][col].getNextReflected(), self.stationaryMirror.m1)):
                    rayBetweenStationPower = Ray(stationaryReflection1, self.degrees)
                    self.rayArray2D_1[row][col].getNextReflected().setNextReflected(rayBetweenStationPower)
                    self.rayArray2D_1[row][col].getNextReflected().reflected = True
                    rayBetweenStationPower.setStartingLocation(stationStartingCoord)
                    rayBetweenStationPower.setLength(lengthToPlane(rayBetweenStationPower, powerTowerSurface, [0,1,0,1])[0])
                else:
                    self.rayArray2D_1[row][col].getNextReflected().setLength(length*100)
                self.rayArray2D_1[row][col].radiallyRotate()
                drawRay(self.rayArray2D_1[row][col])


                self.rayArray2D_4[row][col] = Ray(self.sunRotation, self.degrees)
                startingFour = np.array([startingLocation[0] + helioFourOffset, startingLocation[1], startingLocation[2], 1])
                #Set the starting location on the rotated heliostat for every point
                self.rayArray2D_4[row][col].setStartingLocation(startingFour)
                #Set next ray in the linked list to the reflection of the sun's ray off the rotated heliostat
                self.rayArray2D_4[row][col].setNextReflected(Ray(helioReflection, self.degrees))

                #Set the next ray starting location (it is the same as the sun ray)
                self.rayArray2D_4[row][col].getNextReflected().setStartingLocation(startingFour)
                #Set next ray distance of heliostat to stationary mirror
                length, stationStartingCoord = lengthToPlane(self.rayArray2D_4[row][col].getNextReflected(), self.stationaryMirror.getMirrorFour(), self.stationaryMirror.getNormalFour())
                #Set the next ray length
                self.rayArray2D_4[row][col].getNextReflected().setLength(length)
                #Check if the ray actually hits the stationary mirror
                if(stationaryMirrorRefected(self.rayArray2D_4[row][col].getNextReflected(), self.stationaryMirror.m4)):
                    rayBetweenStationPower = Ray(stationaryReflection4, self.degrees)
                    self.rayArray2D_4[row][col].getNextReflected().setNextReflected(rayBetweenStationPower)
                    self.rayArray2D_4[row][col].getNextReflected().reflected = True
                    rayBetweenStationPower.setStartingLocation(stationStartingCoord)
                    rayBetweenStationPower.setLength(lengthToPlane(rayBetweenStationPower, powerTowerSurface, [0,1,0,1])[0])
                else:
                    self.rayArray2D_4[row][col].getNextReflected().setLength(length*100)
                self.rayArray2D_4[row][col].radiallyRotate()
                drawRay(self.rayArray2D_4[row][col])

    def main(self):

        #Create the window that displays the game, and set initial parameters for the screen
        initialization()

        #Keep track of which line to write to
        currentExcelLine = None

        trackingOn = False
        skyOn = False
        rayOn = True

        self.stationaryMirror.radiallyRotate()
        self.rotatedMirror.radiallyRotate()

        while True:
	        #Slow CPU - 30 screens / sec
            #clock.tick(30)
            for event in pygame.event.get():
                #.QUIT == red X
                if event.type == pygame.QUIT:
                    pygame.quit()
                    quit()
		        #if escape is clicked, quit
                if event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE:
                    pygame.quit()
                    quit()

                #Check for visual coordinate rotation inputs
                self.mouseMove(event)
                self.zoomInOut(event)
                self.absoluteRotate(event)
                self.localAxisRotations(event)

                # If p is pressed, switch boolean for trackingOn
                if(event.type == pygame.KEYDOWN and event.key == pygame.K_p and trackingOn == False):
                    print("\nTracking on\n")
                    trackingOn = True
                elif(event.type == pygame.KEYDOWN and event.key == pygame.K_p and trackingOn):
                    trackingOn = False
                #If r is pressed, switch boolean for rayOn
                if(event.type == pygame.KEYDOWN and event.key == pygame.K_r and rayOn == True):
                    rayOn = False
                elif(event.type == pygame.KEYDOWN and event.key == pygame.K_r and rayOn == False):
                    rayOn = True

            #Clear the current drawing
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

            #Display power tower
            powerTower()

            #Tracking will set a new normal for the heliostat
            if(trackingOn):

                # Contains local x,y rotations for heliostat in radians, vector of sun location
                heliostatRotationArray, sunRotation, solarTime = getTrackingCoordinates()
                self.sunRotation = sunRotation
                if self.solarTime is None:
                    self.initialTime = solarTime
                self.solarTime = solarTime

                #Tracking will draw the rotated heliostat in the position outputted by tracking.py
                helio, helioNormal = self.rotatedMirror.tracking(heliostatRotationArray)

                if currentExcelLine is None:
                    currentExcelLine = 1
                #rayTracing will draw the 'light rays'
                if(rayOn):
                    self.rayTracing(helio, helioNormal)

                    currentTime = self.solarTime.hour*60 + self.solarTime.minute
                #Draw the sun
                sun(self.sunRotation)

            else:

                #Draw the sun if it has been initialized
                if not self.sunRotation is None:
                    sun(self.sunRotation)

                #Draw the heliostat
                self.rotatedMirror.rotatedMirror()

                if self.rayArray2D_2[0][0] is not None and rayOn:
                    for i in range(len(self.rayArray2D_2)):
                        for j in range(len(self.rayArray2D_2[i])):
                            drawRay(self.rayArray2D_1[i][j])
                            drawRay(self.rayArray2D_2[i][j])
                            drawRay(self.rayArray2D_3[i][j])
                            drawRay(self.rayArray2D_4[i][j])

            #display the stationary mirrors
            self.stationaryMirror.stationaryMirror()

            #Displays the local axis lines
            Axis(self.yAbsoluteToLocal, self.xAbsoluteToLocal)

            #Updates the drawing
            pygame.display.flip()

            if currentExcelLine is not None:
                currentExcelLine = currentExcelLine + 1

def initialization():

    X = 1200
    Y = 800

    #Setup pygame
    pygame.init()
    #Size in pixel of display
    display = (X, Y)
    #Initialize display
    screen = pygame.display.set_mode(display, DOUBLEBUF|OPENGL, pygame.RESIZABLE)
    pygame.display.set_caption("Ra Solar Heliostat Visualization")
    #Set the perspective of where the used views the generated field
    glLineWidth(2)
    gluPerspective(45,(display[0]/display[1]), 0.1, 220)
    glTranslatef(0, 0, -15)
    glScaled(0.05, 0.05, 0.05)

def main():
    print("=========================================================================================")
    print("")
    print("Axis:     Local X        Local Y        Local Z        Absolute visual X (right),Y (up)")
    print("          Green           Red       Blue (Toward Sky)              White\n")
    print("Note: This visualization only currently works when the sun is above the horizon, which shouldn't effect anything but can have wierd outputs\n")

    Game().main()

if __name__ == '__main__':
    main()
