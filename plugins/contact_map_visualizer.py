#!/usr/bin/env python

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you are
# free to use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# -----------------------------------------------------------------------------------
# This PyMOL Plugin Contact Maps Visualizer is
# Copyright (C) 2012 by Venkatramanan Krishnamani <venks@andrew.cmu.edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
#

###############################################################################
#                        Contact Maps Visualizer 1.0                          #
###############################################################################

import Tkinter
import tkFileDialog
import Image
import pygame as pg
from pymol import cmd, stored
import os, sys, urllib2, zlib
import random

count = 0
colors = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan', 'orange', 'marine', 'chartreuse', 'limon']
colors_value = [(255,0,0), (0,0,255), (0,255,0), (255, 255, 0), (255,0,255), (0, 255, 255), (255, 128, 0), (0, 128,255), (128, 255, 0), (191, 255, 64)]

def __init__(self):
        self.menuBar.addmenuitem('Plugin', 'command','Contact Map Visualizer',label = 'Contact Map Visualizer',command = lambda imLoc=self : getImageandPDBLocations())
    
def getImageandPDBLocations():
        dialog = Tkinter.Tk()
        dialog.withdraw()
        myFormats = [('Portable Network Graphics','*.png'),('JPEG / JFIF','*.jpg')]
        file = tkFileDialog.askopenfile(parent=dialog,mode='rb',filetypes=myFormats, title='Choose the contact map image file')
        print "Opening...", file.name
        if file != None:
            data = file.read()
            print "Received %d bytes from this file." % len(data)

        myFormatsPDB = [('Protein Data Bank','*.pdb'), ('MDL mol','*.mol'), ('PyMol Session File','*.pse')]
        pdbFile = tkFileDialog.askopenfile(parent=dialog,mode='rb',filetypes=myFormatsPDB, title='Choose the corresponding PDB file')
        print "Opening...", pdbFile.name
        cmd.load(pdbFile.name, pdbFile.name)
        preProcessPDB(pdbFile.name)
        loadImageOnScreen(file.name)

def preProcessPDB(fname):
        stored.residues = []
        stored.chains = []
        cmd.iterate('(name ca)', 'stored.residues.append(resi)')
        cmd.iterate('(name ca)', 'stored.chains.append(chain)')
        cmd.show_as("cartoon", '(all)')
        cmd.color("white", '(all)')
        
def loadImageOnScreen(image_file):
        
        #General variables
        global count
        sel = 0
        rand  = int(random.random()*100000)
        factor = 2
        gray = (100,100,100)
        outputname = "%s_%s.png"%("contact-map-selections",rand)
        print "Output File name : outputname"
        
        #Text related arrays nd variables
        text = []
        textRect = []
        textcount = -1
        BLACK = (0, 0, 0)
        WHITE = (255, 255, 255)

        # use an image you have (.bmp  .jpg  .png  .gif)
        im = Image.open(image_file)
        sw = im.size[0] * factor
        sh = im.size[1] * factor
        
        # initialize pygame
        pg.init()
        screen = pg.display.set_mode((sw, sh))
        pg.display.set_caption('Choose the location by clicking')
        image = pg.image.load(image_file).convert()
        image = pg.transform.scale(image, (sw,sh))
        image_rect = image.get_rect()
        
        running = True
        while running:
            event = pg.event.poll()
            keyinput = pg.key.get_pressed()
            
            # exit on corner 'x' click or escape key press
            if keyinput[pg.K_ESCAPE]:
                pg.quit()
            elif event.type == pg.QUIT:
                running = False
                pg.quit()
            elif event.type == pg.MOUSEBUTTONDOWN:
                sel = 1
                textcount += 1
                coor = [event.pos[0]/2,(sh - event.pos[1])/2]
                updateSelections(coor)
                if count >= 10:
                        count = 0
                pg.draw.circle(image, colors_value[count], event.pos, 3, 0)
                
                # set up fonts
                basicFont = pg.font.SysFont("Arial", 12)
                name = "(%s%s, %s%s)"%(stored.residues[coor[0]], stored.chains[coor[0]],stored.residues[coor[1]], stored.chains[coor[1]])
                # set up the text
                text.append(basicFont.render(name, True, WHITE, BLACK))
                textRect.append(text[textcount].get_rect())
                textRect[textcount][0] = event.pos[0] + 5
                textRect[textcount][1] = event.pos[1] + 5
                
                screen.blit(image, image_rect)
                for a in range(textcount+1):
                        screen.blit(text[a], textRect[a])
                pg.display.flip()
                pg.image.save(screen, outputname)
                count += 1
            # update screen
            screen.blit(image, image_rect)
            if sel == 1:
                for a in range(textcount+1):
                        screen.blit(text[a], textRect[a])
            pg.display.flip()

def updateSelections(residues):
        global count
        print "Residue Numbers of Current Selection (screenshot saved)"
        print "x", stored.residues[residues[0]], stored.chains[residues[0]]
        print "y", stored.residues[residues[1]], stored.chains[residues[1]]
        selectionname = "%s%s_%s%s"%(stored.residues[residues[0]], stored.chains[residues[0]], stored.residues[residues[1]], stored.chains[residues[1]])
        
        cmd.select(selectionname, "resid %s and chain %s + resid %s and chain %s"%(stored.residues[residues[0]],stored.chains[residues[0]],stored.residues[residues[1]], stored.chains[residues[1]]))
        
        cmd.set('stick_radius', '0.2')
        cmd.show("sticks", selectionname)
        cmd.set('sphere_scale', '0.25')
        cmd.show("spheres", selectionname)

        if count >= 10:
                count = 0
        cmd.color(colors[count], selectionname)
