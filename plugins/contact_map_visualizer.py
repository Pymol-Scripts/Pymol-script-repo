'''
Contact Map Visualizer 2.1

http://pymolwiki.org/index.php/Contact_Map_Visualizer

Author: Venkatramanan Krishnamani (Version 1.1)
Author: Thomas Holder (Version 2.0)
'''

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

from __future__ import absolute_import
from __future__ import print_function

import os
from pymol import cmd, CmdException

colors = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan', 'orange',
          'marine', 'chartreuse', 'purpleblue', 'violet', 'limon', ]
colors_value = [tuple(int(i * 255) for i in cmd.get_color_tuple(color))
                for color in colors]


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'Contact Map Visualizer',
                             label='Contact Map Visualizer', command=lambda s=self: CMVDialog(s))


def CMVDialog(self):
    import sys
    if sys.version_info[0] < 3:
        import tkFileDialog
        import tkMessageBox
    else:
        import tkinter.filedialog as tkFileDialog
        import tkinter.messagebox as tkMessageBox

    try:
        import pygame as pg
    except ImportError:
        tkMessageBox.showerror('Error', 'This plugin requires the "pygame" module')
        return

    myFormats = [('Portable Network Graphics', '*.png'), ('JPEG / JFIF', '*.jpg')]
    try:
        image_file = tkFileDialog.askopenfilename(parent=self.root,
                                                  filetypes=myFormats, title='Choose the contact map image file')
        if not image_file:
            raise
    except:
        tkMessageBox.showerror('Error', 'No Contact Map!')
        return

    myFormatsPDB = [('Protein Data Bank', '*.pdb'), ('MDL mol', '*.mol'), ('PyMol Session File', '*.pse')]
    try:
        pdb_file = tkFileDialog.askopenfilename(parent=self.root,
                                                filetypes=myFormatsPDB, title='Choose the corresponding PDB file')
        if not pdb_file:
            raise
    except:
        tkMessageBox.showerror('Error', 'No PDB file!')
        return

    name = cmd.get_unused_name('protein')
    cmd.load(pdb_file, name)
    contact_map_visualizer(image_file, name, 1, 0)


def contact_map_visualizer(image_file='', selection='all', screenshots=0, quiet=1):
    '''
DESCRIPTION

    Contact map visualizer

    If no image_file is provided, try to run "contact_map_generator" first.

USAGE

    contact_map_visualizer [ image_file [, selection [, screenshots ]]]

SEE ALSO

    contact_map_generator
    '''
    if not image_file:
        import tempfile
        image_file = tempfile.mktemp(suffix='.png')
        if not quiet:
            print(' Warning: no image_file provided!')
            print(' Will try to generate it for selection "%s"' % (selection))
            print(' Writing image to', image_file)
        contact_map_generator(image_file, selection, quiet=quiet)
    elif image_file.lower().endswith('.xpm'):
        image_file_png = image_file[:-4] + '.png'
        if not quiet:
            print(' Converting image to', image_file_png)
        xpm_convert(image_file, image_file_png)
        image_file = image_file_png

    import threading
    t = threading.Thread(target=_contact_map_visualizer, kwargs=locals())
    t.setDaemon(1)
    t.start()


def _contact_map_visualizer(image_file, selection, screenshots, quiet, **kwargs):
    from datetime import datetime

    try:
        import pygame as pg
    except ImportError:
        print(' Error: This plugin requires the "pygame" module')
        raise CmdException

    screenshots, quiet = int(screenshots), int(quiet)

    cmd.set('sphere_scale', 0.25)

    idx_list = []
    cmd.iterate('(%s) and name CA' % (selection),
                'idx_list.append(((model,index),chain,resi))', space=locals())
    cmd.color("white", selection)

    # General variables
    count = 0
    sel = 0
    ntime = datetime.now()
    BLACK = (0, 0, 0)
    WHITE = (255, 255, 255)

    if screenshots:
        outputname = "%s_selectedPoints_%d-%d-%d_%d%d.png" % (os.path.basename(image_file),
                                                              ntime.day, ntime.month, ntime.year, ntime.hour, ntime.minute)

    # Text related arrays nd variables
    text = []
    textRect = []
    textcount = -1

    # screen height
    try:
        screen_height = cmd.pymol._ext_gui.root.winfo_screenheight()
    except:
        print(' Warning: could not determine screen height')
        screen_height = 700

    # use an image you have (.bmp  .jpg  .png  .gif)
    image = pg.image.load(image_file)
    image_height = image.get_height()

    factor = max((screen_height - 100) / image_height, 1)
    height = image_height * factor
    factor = 1.0 / factor

    if image_height != len(idx_list):
        print(' Warning: Dimension of image and number of CA atoms in selection', end=' ')
        print('differ! (%d vs. %d)' % (image_height, len(idx_list)))

    # initialize pygame
    pg.init()
    screen = pg.display.set_mode((height, height))
    pg.display.set_caption('Choose the location by clicking')
    image = pg.transform.scale(image.convert(), (height, height))
    image_rect = image.get_rect()

    while True:
        event = pg.event.poll()
        keyinput = pg.key.get_pressed()

        # exit on corner 'x' click or escape key press
        if keyinput[pg.K_ESCAPE] or event.type == pg.QUIT:
            break

        if event.type == pg.MOUSEBUTTONDOWN:
            coor = [event.pos[0] * factor, (height - event.pos[1]) * factor]

            try:
                idx1, chain1, resi1 = idx_list[int(coor[0])]
                idx2, chain2, resi2 = idx_list[int(coor[1])]
            except IndexError:
                print(' Error: selection to small')
                continue

            print(' You clicked %s/%s/ %s/%s/' % (chain1, resi1, chain2, resi2))

            sel = 1
            textcount += 1

            selectionname = '%s%s_%s%s' % (resi1, chain1, resi2, chain2)
            cmd.select(selectionname, 'byres (%s`%d %s`%d)' % tuple(idx1 + idx2))

            cmd.show('sticks', selectionname)
            cmd.show('spheres', selectionname)

            cmd.color(colors[count], selectionname)
            cmd.center(selectionname, animate=1)

            pg.draw.circle(image, colors_value[count], event.pos, 3, 0)

            # set up fonts
            basicFont = pg.font.SysFont('Arial', 12)
            name = '(%s%s, %s%s)' % (resi1, chain1, resi2, chain2)

            # set up the text
            text.append(basicFont.render(name, True, WHITE, BLACK))
            textRect.append(text[textcount].get_rect())
            textRect[textcount][0] = event.pos[0] + 5
            textRect[textcount][1] = event.pos[1] + 5

            screen.blit(image, image_rect)
            for a in range(textcount + 1):
                screen.blit(text[a], textRect[a])
            pg.display.flip()

            if screenshots:
                if not quiet:
                    print(' Writing image to', outputname)
                pg.image.save(screen, outputname)

            count += 1
            if count >= len(colors):
                count = 0

        # update screen
        screen.blit(image, image_rect)
        if sel == 1:
            for a in range(textcount + 1):
                screen.blit(text[a], textRect[a])
        pg.display.flip()

    pg.quit()


def contact_map_generator(filename, selection='all', state=-1, quiet=1):
    '''
DESCRIPTION

    Genarate a contact map image with "g_mdmat" (requires gromacs)

USAGE

    contact_map_generator filename [, selection [, state ]]
    '''
    import os
    import shutil
    import tempfile
    import subprocess

    state, quiet = int(state), int(quiet)

    if not os.path.splitext(filename)[-1].lower() in ['.png', '.jpg', '.jpeg']:
        print(' Error: filename must have png or jpg extension')
        raise CmdException

    try:
        tempdir = tempfile.mkdtemp()
        file_f = os.path.join(tempdir, 'f.pdb')
        file_mean = os.path.join(tempdir, 'dm.xpm')

        cmd.save(file_f, selection, state)

        if state == 0:
            file_s = os.path.join(tempdir, 's.pdb')
            cmd.save(file_s, selection)
        else:
            file_s = file_f

        process = subprocess.Popen(['g_mdmat', '-f', file_f, '-s', file_s, '-mean', file_mean],
                                   stdin=subprocess.PIPE)
        print('Protein-H', file=process.stdin)
        process.stdin.close()
        process.wait()

        xpm_convert(file_mean, filename)

    except OSError as e:
        print(e)
        print(' Error: calling external applications failed')
        raise CmdException
    finally:
        shutil.rmtree(tempdir)


def xpm_convert(infile, outfile):
    '''
    Strips comments and repeated spaces from XPM file and saves it as new file.
    '''
    import re
    try:
        from PIL import Image
    except ImportError:
        import Image
    from io import StringIO

    xpm = open(infile).read()
    xpm = re.sub(r'/\*.*?\*/', '', xpm)  # strip comments
    xpm = re.sub(r'  +', ' ', xpm)      # strip multi-spaces
    xpm = re.sub(r'\n\s+', '\n', xpm)   # strip empty lines
    xpm = u'/* XPM */' + xpm

    image = Image.open(StringIO(xpm))
    image.save(outfile)

cmd.extend('contact_map_visualizer', contact_map_visualizer)
cmd.extend('contact_map_generator', contact_map_generator)

# vi:expandtab:smarttab:sw=4
