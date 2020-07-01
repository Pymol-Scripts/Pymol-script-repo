from __future__ import print_function

import sys
if sys.version_info[0] < 3:
    from Tkinter import *
    import tkSimpleDialog
    import tkMessageBox
    import tkColorChooser
else:
    from tkinter import *
    from tkinter import simpledialog as tkSimpleDialog
    import tkinter.messagebox as tkMessageBox
    import tkinter.colorchooser as tkColorChooser

import re
from pymol import stored, cmd
import math
from pymol.cgo import *
from pymol.vfont import plain

from cctbx import sgtbx, uctbx
import numpy as N
from numpy.linalg import *


def __init__(self):
    #MAIN
    self.menuBar.addcascademenu('Plugin','SuperSym')
    #DEFAULT SET BUILD
    self.menuBar.addmenuitem('SuperSym', 'command', 'Default Symmetry Partner Set',
                             label = 'Default Symmetry Partner Set',
                             command = lambda s = self: symDialog(s, 0))
    #UNIT CELL BUILD
    self.menuBar.addmenuitem('SuperSym', 'command', 'Draw Unit Cell',
                             label = 'Draw Unit Cell',
                             command = lambda s = self: cellDialog(s))
    #SYM SUBMENU
    self.menuBar.addcascademenu('SuperSym', 'Build Symmetry Partners')

    self.menuBar.addmenuitem('Build Symmetry Partners', 'command', 'Cell [0,0,0] (default)',
                             label = 'Cell [0,0,0] (default)',
                             command = lambda s = self: symDialog(s, 0))

    self.menuBar.addmenuitem('Build Symmetry Partners', 'command', 'Cell [x,y,z] (custom)',
                             label = 'Cell [x,y,z] (custom)',
                             command = lambda s = self: symDialog(s, 1))

    self.menuBar.addmenuitem('Build Symmetry Partners', 'command', '2x2x2 Block',
                             label = '2x2x2 Block',
                             command = lambda s = self: symDialog(s, 2))

    self.menuBar.addmenuitem('Build Symmetry Partners', 'command', '3x3x3 Block',
                             label = '3x3x3 Block',
                             command = lambda s = self: symDialog(s, 3))

    self.menuBar.addmenuitem('Build Symmetry Partners', 'command', 'By Partner',
                             label = 'By Partner',
                             command = lambda s = self: symDialog(s, 4))
    #COLOR SUBMENU
    self.menuBar.addcascademenu('SuperSym', 'Coloring')

    self.menuBar.addmenuitem('Coloring', 'command', 'Default Rainbow',
                             label = 'Default Rainbow',
                             command = lambda s = self: colorDialog(s, 0))

    self.menuBar.addmenuitem('Coloring', 'command', 'Select color for each operation',
                             label = 'Select color for each operation',
                             command = lambda s = self: colorDialog(s, 1))

    self.menuBar.addmenuitem('Coloring', 'command', 'Select one color for custom set of operations',
                             label = 'Select one color for custom set of operations',
                             command = lambda s = self: colorDialog(s, 2))
    #GRAPHICS SUBMENU
    self.menuBar.addcascademenu('SuperSym', 'Graphics')

    self.menuBar.addmenuitem('Graphics', 'command', 'Lines',
                             label = 'Lines',
                             command = lambda s = self: graphicsDialog(s, 0))

    self.menuBar.addmenuitem('Graphics', 'command', 'Ribbon',
                             label = 'Ribbon',
                             command = lambda s = self: graphicsDialog(s, 1))
    self.menuBar.addmenuitem('Graphics', 'command', 'Cartoon',
                             label = 'Cartoon',
                             command = lambda s = self: graphicsDialog(s, 2))
    self.menuBar.addmenuitem('Graphics', 'command', 'Sphere Surface (best for printing)',
                             label = 'Sphere Surface (best for printing)',
                             command = lambda s = self: graphicsDialog(s, 3))

    self.menuBar.addmenuitem('Graphics', 'command', 'Surface (high load render)',
                             label = 'Surface (high load render)',
                             command = lambda s = self: graphicsDialog(s, 4))
    #SYM AXES SUBMENU
    self.menuBar.addcascademenu('SuperSym', 'Symmetry Axes')

    self.menuBar.addmenuitem('Symmetry Axes', 'command', 'Build Axes',
                             label = 'Build Axes',
                             command = lambda s = self: axesDialog(s))
    #ADD OTHER SYMMETRY AXES OPTION HERE
    self.menuBar.addmenuitem('SuperSym', 'command', 'Move symmetry partners',
                             label = 'Move symmetry partners',
                             command = lambda s = self: cellShiftInfo(s))
    self.menuBar.addmenuitem('SuperSym', 'command', 'About',
                             label = 'About',
                             command = lambda s = self: aboutInfo(s))
    self.menuBar.addmenuitem('SuperSym', 'command', 'Help',
                             label = 'Help',
                             command = lambda s = self: helpInfo(s))
    cmd.cell_shift = cell_shift
    cmd.get_operations = get_operations
    cmd.get_matrix = get_orthogonalization_matrix
    cmd.symset = symset
    cmd.cell_shift_helper = cell_shift_helper
    cmd.set_key("ALT-6", cell_shift_proxyX1)
    cmd.set_key("ALT-4", cell_shift_proxyX2)
    cmd.set_key("ALT-8", cell_shift_proxyY1)
    cmd.set_key("ALT-2", cell_shift_proxyY2)
    cmd.set_key("ALT-5", cell_shift_proxyZ1)
    cmd.set_key("ALT-1", cell_shift_proxyZ2)


'''
symDialog: Dialog generator and command issuer for generating symmetry partners

This function is called by SuperSymMenu when any symmetry partner generating option is
selected. It creates dialog windows and receives user input for symmetry generation parameters.

@app -- identifies the GUI interface to build dialog boxes onto.
@mode -- determines specific treatment of symmetry building command
'''
def symDialog(app, mode):
    prefix = tkSimpleDialog.askstring('Prefix',
    'Enter desired prefix for these partners:', parent=app.root)
    object = tkSimpleDialog.askstring('Object',
    'Enter object to generate partners from:', parent=app.root)
    if mode == 0: #make default symmetry set in cell [0,0,0]
        symset(prefix, object)
    if mode == 1: #make symmetry set in custom cell
        cell = tkSimpleDialog.askstring('Cell',
        'Enter lattice cell coordinates separated by commas (ex:x,y,z):', parent = app.root)
        x,y,z = cell.split(',')
        x,y,z = int(x),int(y),int(z)
        symset(prefix, object, x, y, z)
    if mode == 2: #make 2x2x2 block of symmetry sets
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    symset(prefix, object, i, j, k)
    if mode == 3: #make 3x3x3 block of symmetry sets
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    symset(prefix, object, i, j, k)
    if mode == 4: #select individual partners by operation and cell
        ops = get_operations(object)
        opString = ""
        for i,x in enumerate(ops):
            opString += str(i) + " : " + x + "\n"
        opIndeces = tkSimpleDialog.askstring("Symmetry Operations", opString +
        "Enter numbers of desired operations separated by commas (ex:0,2,9)", parent = app.root)
        opListStrings = opIndeces.split(",")
        opList = []
        for op in opListStrings:
            opList.append(int(op))
        cell = tkSimpleDialog.askstring('Cell',
        'Enter lattice cell coordinates separated by commas (ex:x,y,z):', parent = app.root)
        x,y,z = cell.split(',')
        x,y,z = int(x),int(y),int(z)
        symset(prefix, object, x,y,z, opList)

'''
colorDialog: Dialog generator for coloring commands

This function colors sets of symmetry partners defined by the user in the
dialog which it generates.

@app -- identifies root menu calling this function
@mode -- determines coloring scheme to execute
'''
def colorDialog(app, mode):
    prefix = tkSimpleDialog.askstring('Prefix',
    'Enter the prefix of symmetry partners to color', parent = app.root)
    if mode == 0: #standard rainbow by symmetry operation
        colors = ["red", "orange", "yellow", "green", "blue", "purple",
              "salmon", "grey", "pink", "teal", "brown", "br0", "aquamarine",
              "deepolive", "dirtyviolet", "slate", "br4", "darksalmon", "br7",
              "chocolate", "firebrick", "brightorange"]
        for i in range(20):
            try: #required because PyMOL inappropriately throws an exception
                 #when the cmd.color() function colors no objects
                cmd.color(colors[i], "%s%02d*" % (prefix, i))
            except:
                pass #allows us to move on to next symmetry operator
    if mode == 1: #specify for each symmetry operation
        cmd.iterate_state(1, prefix + "*", "stored.tmpObject = model")
        ops = get_operations(stored.tmpObject)
        opString = ""
        for i,x in enumerate(ops):
            opString += str(i) + " : " + x + "\n"
        opIndeces = tkSimpleDialog.askstring("Symmetry Operations", opString +
        "Enter numbers of desired operations separated by commas (ex:0,2,9) or all", parent = app.root)
        if opIndeces == "all":
            opList = []
            for i in range(len(ops)):
                opList.append(i)
        else:
            opList = opIndeces.split(",")
        opStringList = opString.split("\n")
        for i in opList:
            try:
                cmd.color("white", "%s%02d*" % (prefix, i))
                cmd.center("%s%02d*" % (prefix, i))

            except:
                pass
            tempColor = tkColorChooser.askcolor(title = "Color for " + opStringList[int(i)] + " (currently white)",
                                                parent = app.root)[0]
            rgb = []
            for value in tempColor:
                value = float(value)
                value = value/255
                rgb.append(value)
            cmd.set_color("tempColor", rgb)
            try:
                cmd.color("tempColor", "%s%02d*" % (prefix, i))
            except:
                pass
    if mode == 2: #monochrome for a set of operations
        cmd.iterate_state(1, prefix + "*", "stored.tmpObject = model")
        ops = get_operations(stored.tmpObject)
        opString = ""
        for i,x in enumerate(ops):
            opString += str(i) + " : " + x + "\n"
        opIndeces = tkSimpleDialog.askstring("Symmetry Operations", opString +
        "Enter numbers of desired operations separated by commas (ex:0,2,9) or all", parent = app.root)
        if opIndeces == 'all':
            opList = []
            for i in range(len(ops)):
                opList.append(i)
        else:
            opList = opIndeces.split(",")
        opStringList = opString.split("\n")
        tempColor = tkColorChooser.askcolor(parent = app.root)[0]
        rgb = []
        for value in tempColor:
            value = float(value)
            value = value/255
            rgb.append(value)
        cmd.set_color("tempColor", rgb)
        for i in opList:
            try:
                cmd.color("tempColor", "%s%02d*" % (prefix, i))
            except:
                pass
'''
graphicsDialog: Dialog generator for graphics commands

This function sets visual representations for sets of symmetry partners.

@app -- identifies root menu
@mode -- determines type of representation to show
'''
def graphicsDialog(app, mode):
    prefix = tkSimpleDialog.askstring('Prefix',
    'Enter prefix of symmetry partners to display', parent = app.root)
    cmd.hide("everything", prefix + "*")
    if mode == 0: # show lines
        cmd.show("lines", prefix + "*")
    if mode == 1: # show ribbon
        cmd.show("ribbon", prefix + "*")
    if mode == 2: # show cartoon
        cmd.show("cartoon", prefix + "*")
    if mode == 3: # sphere surface
        objSel = prefix + "*"
        findSurfaceResidues(objSel, 3.5, "surface")
        cmd.set("sphere_scale", 1.8)
        cmd.show("spheres", "surface")
    if mode == 4: # regular surface
        cmd.show("surface", prefix + "*")

'''
cellDialog: dialog proxy for draw_cell

This function generates a unit cell representation
FUTURE IMPLEMENTATIONS: select which lattice coordinates to generate unit cell for

@app -- identifies root menu
'''
def cellDialog(app):
    object = tkSimpleDialog.askstring('Object',
    'Enter object to generate cell for:', parent = app.root)
    if tkMessageBox.askyesno('3D Printing', 'Going to print this model?', parent = app.root):
        draw_cell(object, 3.0)
    else:
        draw_cell(object)

'''
axesDialog: dialog proxy for draw_symops_cctbx

This function generates one set of symmetry axes for a given object
FUTURE IMPLEMENTATIONS: select individual axes to generate, attach to model for 3D printing,
                        generate axes for multiple unit cells

@app -- identifies root menu
'''
def axesDialog(app):
    object = tkSimpleDialog.askstring('Object',
    'Enter object to generate symmetry axes for:', parent = app.root)
    if tkMessageBox.askyesno('3D Printing', 'Going to print this model?', parent = app.root):
        draw_symops(object, 2.0)
    else:
        draw_symops(object)

'''
cellShiftInfo: displays info for using cell_shift hotkeys

@app -- identifies root menu
'''
def cellShiftInfo(app):
    tkMessageBox.showinfo('Cell Shifting',
    "To shift a symmetry partner, simply click to select any part of it (select only one partner at a time). \n\n" +
    "Next, hold ALT and press the numpad key corresponding to the axis direction you\'d like to move. \n\n" +
    "Key assignments:\n" +
    "A (x) axis: down--4, up--6 \n" +
    "B (y) axis: down--2, up--8 \n" +
    "C (z) axis: down--1, up--5", parent = app.root)
    tkMessageBox.showwarning('Caution', 'Only attempt to shift symmetry partners created by SuperSym.'+
                             'Attempting to shift any other object will result in errors.')

def aboutInfo(app):
    tkMessageBox.showinfo('About',
                          'SuperSym \nDeveloped by Stuart Ballard (srballard@wisc.edu)\nDepartment of Biochemistry\n'+
                          'University of Wisconsin-Madison', parent = app.root)
def helpInfo(app):
    tkMessageBox.showinfo('Help',
                          'For documentation see http://pymolwiki.org/index.php/SuperSym', parent = app.root)

'''
symset: generates up to one full set of symmetry partners for a given object in a given lattice position

1. Obtain all essential symmetry information from CCTBX. This includes the space group, unit cell parameters,
and fractional coordinates corresponding to symmetry operations.
2. Generate transformation matrices to translate coordinates from orthogonal to fractional, and back.
3. 
'''
def symset(prefix = "sym", object = -1, x=0,y=0,z=0, opList = []):
    if object == -1:
        object = cmd.get_names()[0]
    cell = [float(x),float(y),float(z)]
    view = cmd.get_view()
    cmd.show("lines", object)
    sgInfo = cmd.get_symmetry(object)
    raw_ops = []
    for s in sgtbx.space_group_info(sgInfo[6]).group():
        raw_ops.append(str(s))
    if len(opList) == 0:
        for i in range(len(raw_ops)):
            opList.append(i)
    opMatrices = []
    vars = ["x","y","z"]
#CREATE 4X4 MATRICES FOR SYMMETRY OPERATORS
    for i,raw_op in enumerate(raw_ops):
        ops = raw_op.split(",")
        matrix = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]]
        for j in range(len(ops)):
            for k in range(len(vars)):
                index = ops[j].find(vars[k])
                if index != -1:
                    if index == 0:
                        matrix[k][j] = 1
                    elif ops[j][index - 1] == "-":
                        matrix[k][j] = -1
                    else:
                        matrix[k][j] = 1
            index = ops[j].find("/")
            if index != -1:
                matrix[3][j] = float(ops[j][index - 1]) / float(ops[j][index + 1])
        opMatrices.append(matrix)
    a,b,c,alpha,beta,gamma = sgInfo[0:6]
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    sb = math.sin(math.radians(beta))
    sg = math.sin(math.radians(gamma))
    fracToOrt = N.array([[a, b * cg, c * cb, 0.0],
                                [0.0, b * sg, c * (ca - cb * cg) / sg, 0.0],
                                [0.0, 0.0, c * sb * math.sqrt(1.0 - ((cb * cg - ca) / (sb * sg))**2), 0.0],
                                [0.0,0.0,0.0,1.0]])
    fracToOrt = fracToOrt.transpose()
    ortToFrac = inv(fracToOrt)
    stored.atoms = []
    cmd.iterate_state(1,object,"stored.atoms.append([x,y,z,1])")
    stored.atoms = N.array(stored.atoms)
    fracCoords = N.dot(stored.atoms,ortToFrac)
    for i in opList:
        try:
            op = opMatrices[i]
        except:
            print("Bad symmetry partner numbers. Try again.")
            quit()
        copy = "%s%02d_%d_%d_%d" % (prefix, i, x, y, z)
        cmd.copy(copy, object)
        newCoordsFrac = N.dot(fracCoords, op)
        stored.newCoords = N.dot(newCoordsFrac, fracToOrt)
        stored.j = 0
        cmd.alter_state(1,copy,"x,y,z = stored.newCoords[stored.j][0], stored.newCoords[stored.j][1], stored.newCoords[stored.j][2]; stored.j = stored.j + 1")
        xSum=ySum=zSum=0.0
        for a,b,c in newCoordsFrac:
            xSum += a
            ySum += b
            zSum += c
        center = N.array([xSum,ySum,zSum])
        center = center/len(stored.newCoords)
        shift = [cell[0]-math.floor(center[0]),
                 cell[1]-math.floor(center[1]),
                 cell[2]-math.floor(center[2])]
        cell_shift(copy,shift[0],shift[1],shift[2],0)
        '''
        #COPIES COORDINATES OF EACH ATOM TO CORRESPONDING ONE IN GIVEN SYMMETRY PARTNER
        #cmd.alter_state(1, copy, "x,y,z = cmd.sym_partner([x,y,z], stored.tmpOp)")
        #MOVES SYMMETRY PARTNER TO PROPER LATTICE COORDINATES AND CORRECTS FOR NATIVE LATTICE POSITION ERROR
        #stored.xSum,stored.ySum,stored.zSum = 0.0,0.0,0.0
        #atoms = cmd.count_atoms(copy)
        #cmd.iterate_state(1, copy, "stored.xSum = stored.xSum + x; stored.ySum = stored.ySum + y; stored.zSum = stored.zSum + z")
        #xMean = stored.xSum / atoms
        #yMean = stored.ySum / atoms
        #zMean = stored.zSum / atoms
        #xError, yError, zError = N.dot(N.array([xMean,yMean,zMean]), stored.ortToFrac)
        #dX,dY,dZ = cell[0]-math.floor(xError), cell[1]-math.floor(yError), cell[2]-math.floor(zError)
        #cell_shift(copy,dX,dY,dZ, 0)
        '''
    cmd.hide("everything", object)
    cmd.set_view(view)

'''
def sym_partner(coords, op):
    fracCoords = N.dot(N.array(coords), stored.ortToFrac)
    op = op.replace("x", "(" + str(fracCoords[0]) + ")")
    op = op.replace("y", "(" + str(fracCoords[1]) + ")")
    op = op.replace("z", "(" + str(fracCoords[2]) + ")")
    op = op.split(",")
    for i in range(3):
        index = op[i].find("/")
        if index != -1:
            if len(op[i]) == index + 2:
                op[i] = op[i][0:index - 1] + str(float(op[i][index - 1]) / float(op[i][index + 1]))
            else:
                op[i] = op[i][0:index - 1] + str(float(op[i][index - 1]) / float(op[i][index + 1])) + op[i][index + 2:]
        op[i] = eval(op[i])
    return N.dot(N.array(op), stored.fracToOrt)
'''


def cell_shift_proxyX1():
    cmd.iterate_state(1, "sele", "stored.tmpObject = model")
    cell_shift(stored.tmpObject, 1,0,0)
def cell_shift_proxyX2():
    cmd.iterate_state(1, "sele", "stored.tmpObject = model")
    cell_shift(stored.tmpObject, -1,0,0)
def cell_shift_proxyY1():
    cmd.iterate_state(1, "sele", "stored.tmpObject = model")
    cell_shift(stored.tmpObject, 0,1,0)
def cell_shift_proxyY2():
    cmd.iterate_state(1, "sele", "stored.tmpObject = model")
    cell_shift(stored.tmpObject, 0,-1,0)
def cell_shift_proxyZ1():
    cmd.iterate_state(1, "sele", "stored.tmpObject = model")
    cell_shift(stored.tmpObject, 0,0,1)
def cell_shift_proxyZ2():
    cmd.iterate_state(1, "sele", "stored.tmpObject = model")
    cell_shift(stored.tmpObject, 0,0,-1)

def cell_shift(object, dX, dY, dZ, rename = 1):
    if rename:
        oldName = object.split("_")
        oldPre = oldName[0]
        oldX = int(oldName[1])
        oldY = int(oldName[2])
        oldZ = int(oldName[3])
        newX = "_" + str(int(dX) + oldX)
        newY = "_" + str(int(dY) + oldY)
        newZ = "_" + str(int(dZ) + oldZ)
        newName = oldPre + newX + newY + newZ
        #if cmd.get_names().find(newName) != -1:
        #    print "Symmetry partner already exists in destination position!"
        #    quit()
        cmd.set_name(object, newName)
        object = newName
    stored.shift = [float(dX),float(dY),float(dZ)]
    stored.sgInfo = cmd.get_symmetry(object)
    a,b,c,alpha,beta,gamma = stored.sgInfo[0:6]
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    sb = math.sin(math.radians(beta))
    sg = math.sin(math.radians(gamma))
    stored.fracToOrt = N.array([[a, b * cg, c * cb],
                                [0.0, b * sg, c * (ca - cb * cg) / sg],
                                [0.0, 0.0, c * sb * math.sqrt(1.0 - ((cb * cg - ca) / (sb * sg))**2)]])
    stored.fracToOrt = stored.fracToOrt.transpose()
    stored.ortToFrac = inv(stored.fracToOrt)
    cmd.cell_shift_helper = cell_shift_helper
    cmd.alter_state(1, object, "x,y,z = cmd.cell_shift_helper([x,y,z],stored.shift)")

def cell_shift_helper(coords, shift):
     fracCoords = N.dot(N.array(coords), stored.ortToFrac)
     for i in range(3):
         fracCoords[i] += shift[i]
     coords = N.dot(N.array(fracCoords), stored.fracToOrt)
     return coords[0], coords[1], coords[2]

def get_operations(object):
    raw_ops = []
    sgInfo = cmd.get_symmetry(object)
    for s in sgtbx.space_group_info(sgInfo[6]).group():
        raw_ops.append(str(s))
    return raw_ops

def get_orthogonalization_matrix(object, quiet = 0):
    a,b,c,alpha,beta,gamma = cmd.get_symmetry(object)[0:6]
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    sb = math.sin(math.radians(beta))
    sg = math.sin(math.radians(gamma))
    fracToOrt = N.array([[a, b * cg, c * cb],
                                [0.0, b * sg, c * (ca - cb * cg) / sg],
                                [0.0, 0.0, c * sb * math.sqrt(1.0 - ((cb * cg - ca) / (sb * sg))**2)]])
    if not quiet:
        print(fracToOrt)
        print(inv(fracToOrt))
    return fracToOrt

# -*- coding: utf-8 -*-
def findSurfaceResidues(objSel="(all)", cutoff=2.5, selName = 0):
 """
 findSurfaceResidues
  finds those residues on the surface of a protein
  that have at least 'cutoff' exposed A**2 surface area.

 PARAMS
  objSel (string)
   the object or selection in which to find
   exposed residues
   DEFAULT: (all)

  cutoff (float)
   your cutoff of what is exposed or not.
   DEFAULT: 2.5 Ang**2

  asSel (boolean)
   make a selection out of the residues found

 RETURNS
  (list: (chain, resv ) )
   A Python list of residue numbers corresponding
   to those residues w/more exposure than the cutoff.

 """
 tmpObj="__tmp"
 cmd.create( tmpObj, objSel + " and polymer");
 cmd.set("dot_solvent");
 cmd.get_area(selection=tmpObj, load_b=1)

 # threshold on what one considers an "exposed" atom (in A**2):
 cmd.remove( tmpObj + " and b < " + str(cutoff) )

 stored.tmp_dict = {}
 cmd.iterate(tmpObj, "stored.tmp_dict[(chain,resv)]=1")
 exposed = list(stored.tmp_dict.keys())
 exposed.sort()

 cmd.select(selName, objSel + " in " + tmpObj )
 cmd.delete(tmpObj)

 return exposed

#CELL DRAWING

def set_to_zero(a):
  if abs(a) < 1e-10:
    a=0
  return a

def draw_cell(obj,radius=1.0,mode=0):
  """
  From pymol issue the "run draw_cell.py" command to load the script,
  then issue the "draw_cell(object,<optional radius>)" command
  to actually run it and create the cgo object showing the unit cell
  border for the space group specified by molecular object 'object'.

  e.g. load 1avv.pdb
       run draw_cell.py
       draw_cell 1avv 0.5   (or draw_cell('1avv',.5))

  see also help(draw_cell_param) to draw the cell border for
  user-defined cell dimensions (i.e. not loaded from a pdb file)

  See also "help(draw_cell_param) to draw the cell border by
  specifying the unit cell parameters directly (i.e. not loaded from
  a pdb file).
  """
  radius=float(radius)
  cell_info=cmd.get_symmetry(obj)
  draw_cell_param(cell_info[0:6],radius,mode)

def draw_cell_param(cell_param_list,radius=1.0,mode=0):
  """
  If you wish to draw the unit cell border for any cell without the
  need to load a pdb file, then do this:

  e.g. run draw_cell.py
       draw_cell_param((45.2,45.2,70.8,90.,90.,120.),0.5)

  to generate the cell border for this trigonal space group "p 31 2 1"
  with a radius of 0.5A.  Labels for the origin, and A, B and C axes
  will appear as well.  The perimeter of the cell is colored with the
  RGB components corresponding to the A,B,C components.
  """

  U=uctbx.unit_cell((cell_param_list))

  vert_000 = list(map(set_to_zero,U.orthogonalize((0.,0.,0))))
  vert_100 = list(map(set_to_zero,U.orthogonalize((1.,0.,0))))
  vert_010 = list(map(set_to_zero,U.orthogonalize((0.,1.,0))))
  vert_001 = list(map(set_to_zero,U.orthogonalize((0.,0.,1))))
  vert_110 = list(map(set_to_zero,U.orthogonalize((1.,1.,0))))
  vert_011 = list(map(set_to_zero,U.orthogonalize((0.,1.,1))))
  vert_101 = list(map(set_to_zero,U.orthogonalize((1.,0.,1))))
  vert_111 = list(map(set_to_zero,U.orthogonalize((1.,1.,1))))

#  vert_000 = map(None,U.orthogonalize((0.,0.,0)))
#  vert_100 = map(None,U.orthogonalize((1.,0.,0)))
#  vert_010 = map(None,U.orthogonalize((0.,1.,0)))
#  vert_001 = map(None,U.orthogonalize((0.,0.,1)))
#  vert_110 = map(None,U.orthogonalize((1.,1.,0)))
#  vert_011 = map(None,U.orthogonalize((0.,1.,1)))
#  vert_101 = map(None,U.orthogonalize((1.,0.,1)))
#  vert_111 = map(None,U.orthogonalize((1.,1.,1)))

  #print vert_000

  #CYLINDER = ['CYLINDER']
  #radius = [0.2]
  #print radius
  cell = []
  cell.append(CYLINDER)
  cell.extend(vert_000 + vert_100 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_000 + vert_010 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_000 + vert_001 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_100 + vert_110 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_100 + vert_101 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_010 + vert_110 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_010 + vert_011 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_001 + vert_101 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_001 + vert_011 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_110 + vert_111 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_101 + vert_111 + [radius] + [1,1,1] + [1,1,1])
  cell.append(CYLINDER)
  cell.extend(vert_011 + vert_111 + [radius] + [1,1,1] + [1,1,1])
  cell.append(SPHERE)
  cell.extend(vert_000 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_001 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_010 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_011 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_100 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_101 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_110 + [radius])
  cell.append(SPHERE)
  cell.extend(vert_111 + [radius])

  cmd.load_cgo(cell,"cell")
  #return cell

  if mode == 1:
    text = [COLOR, 1.0, 0.0, 1.0,]

  #wire_text(text,plain,[-5.,-5.,-1],'Origin',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
  #wire_text(text,plain,map(None,U.orthogonalize((1.05,0.0,0.0))),'A',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
  #wire_text(text,plain,map(None,U.orthogonalize((0.0,1.05,0.0))),'B',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
  #wire_text(text,plain,map(None,U.orthogonalize((0.0,0.0,1.05))),'C',[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])

    cyl_text(text,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[1.0,0.0,1.0])
    cyl_text(text,plain,list(U.orthogonalize((1.05,0.0,0.0))),'A',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[1.0,0.0,0.0])
    cyl_text(text,plain,list(U.orthogonalize((0.0,1.05,0.0))),'B',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[0.0,1.0,0.0])
    cyl_text(text,plain,list(U.orthogonalize((0.0,0.0,1.05))),'C',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]],color=[0.0,0.0,1.0])

    cmd.load_cgo(text,'text')


#AXES DRAWING
#! /usr/bin/env python
# Copyright (c) 2004 Robert L. Campbell

#import math

def set_to_zero(a):
  if abs(a) < 1e-10:
    a=0
  return a

def draw_symbol(start,end,symb,color,radius=0.2):
  degtorad = N.pi/180.
  costhirty = N.cos(30.0*degtorad)
  sinthirty = N.sin(30.0*degtorad)
  symb_obj = []

  if symb in ('2', '2^1'):
    pass

  elif symb in ('3', '3^1', '3^2'):
    symb_obj = [ BEGIN, TRIANGLES, COLOR ] + color
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([radius, 0, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([-radius*sinthirty, radius*costhirty, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([-radius*sinthirty, -radius*costhirty, 0]))[0].tolist())

    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([radius, 0, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([-radius*sinthirty, radius*costhirty, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([-radius*sinthirty, -radius*costhirty, 0]))[0].tolist())
    symb_obj.append(END)

  elif symb in ('4', '4^1', '4^2', '4^3'):
    symb_obj = [ BEGIN, TRIANGLES, COLOR ] + color
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([radius, radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([-radius, radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([-radius, -radius, 0]))[0].tolist())

    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([radius, radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([radius, -radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([start]) + N.array([-radius, -radius, 0]))[0].tolist())

    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([radius, radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([-radius, radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([-radius, -radius, 0]))[0].tolist())

    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([radius, radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([radius, -radius, 0]))[0].tolist())
    symb_obj.append(VERTEX)
    symb_obj.extend((N.array([end]) + N.array([-radius, -radius, 0]))[0].tolist())
    symb_obj.append(END)

  elif symb in ('6', '6^1', '6^2', '6^3', '6^4', '6^5'):
    # hexagons still need to be created :)
    pass

  return symb_obj

def draw_symops(obj,radius=0.2,extension=0):
  """
  From pymol issue the "run draw_symops_cctbx.py" command to load the script,
  then issue the "draw_symops(object,<optional radius>,<optional extension>)" command
  to actually run it and create the cgo object.

  e.g. load 1avv.pdb
       run draw_symops_cctbx.py
       draw_symops 1avv, 0.5, .2
         or draw_symops('1avv',.5,.2)
         or draw_symops 1avv, radius=.5, extension=.2

  The different axis types appear as different objects on the PyMOL menu so they can be turned
  on and off individually.

  See also help(draw_symops_param) to draw operators by specifying the space group
  and cell dimensions directly (i.e. not loaded from a pdb file)

  The 'extension' parameter is a fractional increase in the length of each symmetry
  operator axis drawn.  i.e. a value of 0 is the default and a value of .2 increases
  the length by 20% at each end
  """
  radius=float(radius)
  extension=float(extension)
  cell_info=cmd.get_symmetry(obj)
  draw_symops_param(cell_info[0:6],cell_info[6],radius,extension)

def draw_symops_param(cell_param_list,sg,radius=0.2,extension=0):
  """
  If you wish to draw the symmetry operators for any cell without the need to load a
  pdb file, then do this:

  e.g. run draw_symops_cctbx.py
       draw_symops_param((45.2,45.2,70.8,90.,90.,120.),'p3121',0.5,0.1)

  to generate the symmetry operators for this trigonal space group "p 31 2 1"
  of radius .5 with 10% added as an extension at each end.
  """
  radius=float(radius)
  extension=float(extension)

  U=uctbx.unit_cell((cell_param_list))

#rotation axes
#    "2" "yellow",
#    "3" "orange",
#    "4" "mauve",
#    "6" "purple",

#screw axes (all sub_1 axes are green)
#    "21" "green",
#    "31" "green",
#    "32" "lime",
#    "41" "green",
#    "42" "cyan",
#    "43" "iceblue",
#    "61" "green",
#    "62" "silver",
#    "63" "cyan",
#    "64" "iceblue",
#    "65" "blue",

  color = {
    "2" : [1.0, 1.0, 0.0],
    "3" : [1.0, 0.5, 0.0],
    "4" : [1.0, 0.5, 1.0],
    "6" : [1.0, 0.0, 1.0],
    "2^1" : [0.0, 1.0, 0.0],
    "3^1" : [0.0, 1.0, 0.0],
    "3^2" : [0.5, 1.0, 0.5],
    "4^1" : [0.0, 1.0, 0.0],
    "4^2" : [0.0, 1.0, 1.0],
    "4^3" : [0.5, 0.5, 1.0],
    "6^1" : [0.0, 1.0, 0.0],
    "6^2" : [0.8, 0.8, 0.8],
    "6^3" : [0.0, 1.0, 1.0],
    "6^4" : [0.5, 0.5, 1.0],
    "6^5" : [0.0, 0.0, 1.0],
    }

  sg = sg.upper()
  symop_axes = get_all_axes(sg,extension=extension)

  #CYLINDER = 'CYLINDER'
  ax_obj = {}
  #vert_obj = []

  #debug_out = open('debug.log','w')

  if symop_axes:
    for ax in symop_axes:
      #print ax
      start = list(map(set_to_zero,U.orthogonalize(list(ax['start']))))
      end = list(map(set_to_zero,U.orthogonalize(list(ax['end']))))
###############################################################################
# Tried rounding off start and end values in order to understand why axes go
# missing in the drawing, but seem to be present in the cgo.  Doesn't help!
# e.g. for space group 'p23' one of the 3-fold rotations is missing (0,0,0 -> x,-x,x)
# changing one cell axis to something ever so slightly different recovers the axis
# e.g. set cell to be (30.00001,30.,30.,90.,90.,90) and it works!
#    start = map(lambda x: round(x,3),U.orthogonalize(ax['start']))
#    end = map(lambda x: round(x,3),U.orthogonalize(ax['end']))
###############################################################################
      symb_ax = ax['symb']
      color_ax = color[symb_ax]

      #print "axis: ",symb_ax, start, end
      if symb_ax in ax_obj:
        ax_obj[symb_ax].append(CYLINDER)
      else:
        ax_obj[symb_ax] = [CYLINDER]

      ax_obj[symb_ax].extend(start + end + [radius])
      ax_obj[symb_ax].extend(color[symb_ax] + color[symb_ax])
      ax_obj[symb_ax].extend(draw_symbol(start,end,symb_ax,color[symb_ax],radius*6.))

#  #######################################################################################
#    # Debugging output to try to understand why some axes go missing in the drawing.
#    # They don't appear to be missing from the cgo object, though!
#    for xxx in ax_obj[symb_ax]:
#      if xxx == 9.0:
#        #print "\n\n",xxx
#        xxx = "\n\n" + str(xxx) + " "
#        debug_out.write(xxx)
#      else:
#        #print xxx
#        #xxx = "\n" + str(xxx) + " "
#        xxx = str(xxx) + " "
#        debug_out.write(xxx)
#      #print ax_obj[symb_ax]
#  debug_out.write("\n\n")
#  big_string = str(ax_obj)
#  debug_out.write(big_string)
#  # End of debugging output
#  #######################################################################################

  else:
    print("\nNo symmetry axes found for this space group: %s\n" % sg)

  for key,val in ax_obj.items():
    name=sg + "_" + key
    cmd.load_cgo(val,name)
    #debug_out.write("\n\n" + key + "\n" + str(val))
  #return ax_obj

#cmd.extend("draw_symops_param",draw_symops_param)

#! /usr/bin/env python
# List all axes in the unit cell.

# usage:
#   python all_axes.py     - show axes for the 230 reference settings.
#   python all_axes.py P2  - show axes for (e.g.) space group P2

# RWGK = Ralf W. Grosse-Kunstleve
# RWGK Some further refinement is required:
# RWGK   - List only the axes of highest order (e.g. only 4, not 4 and 2).
# RWGK   - List only the axes with the smallest intrinsic component
# RWGK     (e.g. list only 3(1), not both 3(1) and 3(2)).
# RWGK See also: comment regarding shift_range below.

def list_plus(lhs, rhs):
  return [l + r for l, r in zip(lhs, rhs)]

###def fract_2_dec(fraction):
###  list = fraction.split('/')
###  if len(list) == 2 and list[1] != 0:
###    decimal = float(list[0])/float(list[1])
###  else:
###    decimal = float(fraction)
###  return decimal

def rlc_RTMxAnalysis(M):
  r_info = sgtbx.rot_mx_info(M.r())
  t_info = sgtbx.translation_part_info(M)
  t_intrinsic = t_info.intrinsic_part().mod_positive().as_double()
  t_shift = t_info.origin_shift().mod_positive().as_double()

  #End = list_plus(Start + map(None,r_info.ev()))
####debug
###  trans = 0
###  length = 0
####debug

  #if r_info.type() == 1:
  if r_info.type() < 2:
    #(rt, start, end) = ('1',(0,0,0),(0,0,0))
    return None
  #elif r_info.type() == -1:
  #  (rt, start, end) = (str(r_info.type()),t_shift,())
  elif abs(r_info.type()) == 2:
    trans = sum(t_intrinsic)
    if trans == 0:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      (rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r_info.ev())))
    else:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      (rt, start, end) = (str(r_info.type())+"^1",t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type())+"^1",t_shift,tuple(list_plus(t_shift,r_info.ev())))
  elif r_info.type() == 3:
    if r_info.sense() >= 0 :
      # ignore opposite sense of rotation axes since they superimpose
      trans = N.sqrt(sum((map(lambda x,y:(y-x)*(y-x),(0,0,0),t_intrinsic))))
#      trans = N.sqrt(t_intrinsic[0]**2 + t_intrinsic[1]**2 + t_intrinsic[2]**2)
      if trans == 0:
        maxr = max([abs(x) for x in r_info.ev()])
        r = [float(x)/maxr for x in r_info.ev()]
# fudge to make sure that PyMOL actually draws the axis (move it slightly off [1,-1,1]) !!!
        r[0] = r[0]*1.000001
        (rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r)))
        #(rt, start, end) = (str(r_info.type()),t_shift, tuple(list_plus(t_shift,r_info.ev())))
      else:
        maxr = max([abs(x) for x in r_info.ev()])
        r = [float(x)/maxr for x in r_info.ev()]
        #(rt, start, end) = (str(r_info.type())+ "^" + subscript ,t_shift,tuple(list_plus(t_shift,r)))
        (start, end) = (t_shift,tuple(list_plus(t_shift,r)))
        length = N.sqrt(sum((map(lambda x,y:(y-x)*(y-x),start, end))))

#  r_info.sense() for 3^1 and 3^2 seems always to be "1" ???
#        if r_info.sense() < 0:
#          subscript = str(1-r_info.sense())
#        else:
#          subscript = str(r_info.sense())

# use ratio of trans to length to get the correct axis symbol:
# fudged the value to get the right numbers. (using length/2., rather than length/3.)
        if trans < length*0.5 :
          subscript = '1'
        else:
          subscript = '2'

        rt = str(r_info.type())+ "^" + subscript
        #(rt, start, end) = (str(r_info.type()) + "^" + subscript,t_shift, tuple(list_plus(t_shift,r_info.ev())))
###        print "Type, sense, Start, End, length, trans", rt, r_info.sense(), start, end, length, trans
#        print "type: %s, sense: %s, trans: %s, length: %s," % (r_info.type(), r_info.sense(), trans, length)
#        print "(rt, start, end)", (rt,start,end)
    else:
      return None
  #return (r_info.type(),r_info.ev(), t_intrinsic, t_shift)
  elif r_info.sense() > 0:
    # ignore opposite sense of rotation axes since they superimpose
    trans = sum(t_intrinsic)
    if trans == 0:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      (rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type()),t_shift, tuple(list_plus(t_shift,r_info.ev())))
    else:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      subscript =  str(int(trans*r_info.type()+.5))  # add 0.5 to fix rounding errors
      (rt, start, end) = (str(r_info.type())+ "^" + subscript ,t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type()) + "^" + subscript,t_shift, tuple(list_plus(t_shift,r_info.ev())))
  #return (r_info.type(),r_info.ev(), t_intrinsic, t_shift)
  else:
    return None
#  print "type: %s, sense: %s, trans: %s, length: %s," % (r_info.type(), r_info.sense(), trans, length),
#  print "(rt, start, end)", (rt,start,end)
  return (rt, start, end)

def get_all_axes(space_group_symbol=None, space_group_info=None, extension=0):
  assert space_group_symbol is None or space_group_info is None
  shift_range = 1 # RWGK Works for the 230 reference settings; it is not
          # RWGK clear to me (rwgk) what value is needed in general.
  if space_group_symbol is not None:
    space_group_info = sgtbx.space_group_info(symbol=space_group_symbol)
  #space_group_info.show_summary()

  axes_dict = {}
  for smx in space_group_info.group():
    r = smx.r()
    t = smx.t()
    shift = [0,0,0]
    for shift[0] in range(-shift_range,shift_range+1):
      for shift[1] in range(-shift_range,shift_range+1):
        for shift[2] in range(-shift_range,shift_range+1):
          ts = t.plus(sgtbx.tr_vec(shift, 1)).new_denominator(t.den())
          m = sgtbx.rt_mx(r, ts)
          #print m
          rtmxanal = rlc_RTMxAnalysis(m)
          #print r, t, shift, ts, m
          if rtmxanal:
            #print rtmxanal
            axes_dict[rtmxanal] = 0
  axes_list = list(axes_dict.keys())
  axes_list.sort()

  # reject nonenantiomorphic space groups
  if len(axes_list) > 0 and not re.compile("[A-z]").search(space_group_symbol[1:]):
    try:
      sgtbx.space_group_info(space_group_symbol).show_summary(),
      #print len(axes_list), space_group_symbol
    except:
      print(space_group, space_group_symbol)
      print()
      sys.exit(1)
    axes = []
    for a in axes_list:
      if len(a) == 3 and len(a[1]) == 3 and len(a[2]) == 3:
        tmp_dict = {}
        print("%4s %7.4f %7.4f %7.4f    %7.4f %7.4f %7.4f " % (a[0],a[1][0],a[1][1],a[1][2],a[2][0],a[2][1],a[2][2]))
        tmp_dict['symb'] = a[0]
        start_array = N.asarray(a[1])
        end_array = N.asarray(a[2])
        start_vec = start_array - (end_array - start_array)*extension
        end_vec = end_array + (end_array - start_array)*extension
        tmp_dict['start'] = start_vec
        tmp_dict['end'] = end_vec
#rlc#        tmp_dict['start'] = a[1]
#rlc#        tmp_dict['end'] = a[2]
        axes.append(tmp_dict)
      else:
        print(a)
  else:
    return None

  return axes

if __name__ == "__main__":
  import sys
  if len(sys.argv) == 1:
    for i in range(230):
      get_all_axes(i + 1)
  else:
    for symbol in sys.argv[1:]:
      get_all_axes(symbol)
