# -*- coding: utf-8 -*-
from pymol.cgo import *
from pymol import cmd
import itertools
from random import randint

#############################################################################
#
# drawgridbox.py -- Draw a grid box around a selection
#
# AUTHOR: Cunliang Geng
# DATE  : 19/5/2016
#
# Acknowledgement:
# This scirpt was written based on the DrawBoundingBox by Jason Vertrees
#
#############################################################################
def drawgridbox(selection="(all)", nx=10, ny=10, nz=10, padding=0.0, lw=2.0, r=1.0, g=1.0, b=1.0):
    """
    DESCRIPTION
        Given selection, draw a grid box around it.

    USAGE:
        drawgridbox [selection, [nx, [ny, [nz, [padding, [lw, [r, [g, b]]]]]]]]

    PARAMETERS:
        selection,    the selection to enboxen
                      defaults to (all)

        nx,           number of grids on axis X
                      defaults to 10

        ny,           number of grids on axis Y
                      defaults to 10

        nz,           number of grids on axis Z
                      defaults to 10

        padding,      defaults to 0

        lw,           line width
                      defaults to 2.0

        r,            red color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        g,            green color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        b,            blue color component, valid range is [0.0, 1.0]
                      defaults to 1.0

    RETURNS
        string, the name of the CGO box

    NOTES
        * This function creates a randomly named CGO grid box. The user can
        specify the number of grids on X/Y/Z axis, the width of the lines,
        the padding and also the color.
    """

    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

    print("Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))

    minX = minX - float(padding)
    minY = minY - float(padding)
    minZ = minZ - float(padding)
    maxX = maxX + float(padding)
    maxY = maxY + float(padding)
    maxZ = maxZ + float(padding)
    nX=int(nx)
    nY=int(ny)
    nZ=int(nz)
    dX = (maxX-minX)/nX
    dY = (maxY-minY)/nY
    dZ = (maxZ-minZ)/nZ

    if padding != 0:
        print("Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))
    gridbox = [
        LINEWIDTH, float(lw),
        BEGIN, LINES,
        COLOR, float(r), float(g), float(b),
        ]
    for i in range(nX):
        for j in range(nY):
            for k in range(nZ):
                dots= [
                    VERTEX, minX+i*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                ]
                gridbox += dots
    gridbox.append(END)

    boxName = "gridbox_" + str(randint(0,10000))
    while boxName in cmd.get_names():
        boxName = "gridbox_" + str(randint(0,10000))

    cmd.load_cgo(gridbox,boxName)
    return boxName

cmd.extend ("drawgridbox", drawgridbox)