## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#!/usr/bin/env python
#
# $Id: pixelMap2D.py,v 1.4 2007/10/08 18:16:25 rhuey Exp $

from DejaVu.colorTool import TkColor, RedWhiteRamp
from DejaVu.colorMap import ColorMap
import Image, numpy.oldnumeric as Numeric


class PixelMap2D:
    """
    PixelMap2D:
        Used to display 2D representations of 2D Numeric arrays of data, colored by magnitude.
        Its cmap, a DejaVu.colorMap, provides a mapping from each data
        value to an entry in a 'ramp' which is a list of lists of 3 floats between 0-1. 
        (By default, a RedWhiteRamp is used. Each entry defines a color for that value.)  
        The property range is set by default from the 'mini' and 'maxi' values of the
        results array but it can be set by specifying values for 'mini' and 'maxi'. In
        the latter case, the values of the result array are clamped to be
        greater than 'mini' and less than 'maxi'. A copy of the original data is used.

        Its method to display data "show_image" -> [x,y] produces an array of pixels 
        colored by magnitude of corresponding x,y result and is implemented using 
        PIL.Image show_image.
        "On Unix platforms, this saves the image to a temporary PPM file, and
        calls the xv utility.  On Windows, it saves the image to a temporary BMP file, 
        and uses the standard BMP display utility to show it.""  
        (F.Lundh, http://effbot.org/imagingbook/image.htm)

        Its method to display the corresponding color map legend "show_legend"
        produces an array of pixels, 1-by-nsteps, showing the range of the
        colormap.  This legend is displayed using xv on Unix platforms
    """

    def __init__(self, results, master=None, mini=None, maxi=None, 
                    cmap=None, cmap_name='cmap', ramp=None, npixels=10):
        assert isinstance(results, type(Numeric.array(range(4))))
        assert len(results.shape)==2
        # a copy of the input array is used
        self.results = Numeric.array(results[:])
        self.results = results
        self.master = master
        #for a Numeric array,
        self.x_number = results.shape[0]   #number of rows 
        self.y_number = results.shape[1]   #number of columns
        self.results = results
        self.npixels = npixels
        if mini is None:
            mini=min(results.ravel())
        else:
            for x in xrange(self.x_number):
                for y in xrange(self.y_number):
                    if self.results[x][y]<mini:
                        self.results[x][y]=mini
        self.mini = mini
        if maxi is None:
            maxi=max(results.ravel())
        else:
            for x in xrange(self.x_number):
                for y in xrange(self.y_number):
                    if self.results[x][y]>maxi:
                        self.results[x][y]=maxi
        self.maxi = maxi
        self.cmap = cmap 
        if cmap is None:
            if ramp is None:
                ramp = RedWhiteRamp()
            self.cmap = ColorMap(cmap_name, mini=mini, maxi=maxi, ramp=ramp)

    def show_image(self):
        #'RGB' requires a tuple (r,g,b) where r,g,b are ints in range 0-255
        #'RGBA' requires a tuple (r,g,b,a) where r,g,b,a are ints in range 0-255
        values = self.cmap.Map(self.results.ravel())
        self.image = Image.fromstring('RGBA', (self.x_number,self.y_number), 
                                (values*255).astype('B').tostring())
        self.image.show()


    def show_legend(self, steps=20, horizontal=True):
        delta = (self.maxi-self.mini)/(steps-1)
        vals = self.vals= []
        for i in range(steps-1):
            vals.append(self.mini+i*delta)
        self.vals.append(self.maxi)
        values = self.cmap.Map(vals)
        if horizontal:
            self.legend = Image.fromstring('RGBA', (steps, 1), (values*255).astype('B').tostring())
        else:
            self.legend = Image.fromstring('RGBA', (1,steps), (values*255).astype('B').tostring())
        self.legend.show()


      
if __name__=='__main__':
    import numpy.oldnumeric as Numeric
    from time import sleep
    r = Numeric.array(range(10000))
    r.shape = (100,100)

    rw_w = PixelMap2D(r[:])
    rw_w.show_image()
    print "rw_w shown"
    sleep(3)


    rw_half_w = PixelMap2D(Numeric.array(r[:])[:], maxi=5000.)
    rw_half_w.show_image()
    print "rw_half_w shown"
    sleep(3)

    from DejaVu.colorTool import RedWhiteBlueRamp
    rwb_w = PixelMap2D(r[:], ramp=RedWhiteBlueRamp())
    rwb_w.show_image()
    print "rwb_w shown"
    sleep(3)

    s_r = Numeric.array(Numeric.swapaxes(r[:],0,1))
    s_r.shape = (100,100)
    s_r_w = PixelMap2D(s_r, ramp=RedWhiteBlueRamp())
    s_r_w.show_image()
    print "s_r_w shown"
    sleep(3)

    
    rwb_half_w = PixelMap2D(Numeric.array(r[:])[:], maxi=5000., ramp=RedWhiteBlueRamp())
    rwb_half_w.show_image()
    print "rwb_half_w shown"
    sleep(3)


## this example depends on finding data in working_directory_file 'ligand_results.py'
## stored as a dictionary 'ligand_results'
##stems = ['1a30', '1a94', '1a9m', '1aaq', '1ajv', '1ajx', '1b11', '1b6j', '1b6k',
##    '1b6l', '1b6m', '1b6n', '1b6o', '1b6p', '1bdl', '1bdq', '1bdr', '1bv7',
##    '1bv9', '1bwa', '1bwb', '1d4k', '1d4l', '1dmp', '1fiv', '1fmb', '1g2k',
##    '1g35', '1gnm', '1gnn', '1gno', '1hbv', '1hef', '1heg', '1hih', '1hii',
##    '1hiv', '1hos', '1hpo', '1hps', '1hpv', '1hpx', '1hsg', '1hsh', '1hte',
##    '1htf', '1htg', '1hvh', '1hvi', '1hvj', '1hvk', '1hvl', '1hvr', '1hvs',
##    '1hwr', '1hxw', '1iiq', '1izh', '1izi', '1jld', '1k6c', '1k6p', '1k6t',
##    '1k6v', '1kzk', '1mes', '1met', '1meu', '1mtr', '1mui', '1ody', '1pro',
##    '1qbr', '1qbs', '1qbt', '1qbu', '1sbg', '1tcw', '1tcx', '1vij', '1vik',
##    '2bpv', '2bpx', '2fmb', '2mip', '2upj', '3aid', '3tlh', '4fiv', '4hvp',
##    '4phv', '5fiv', '5hvp', '5upj', '6fiv', '6upj', '7hvp', '7upj', '8hvp', '9hvp']
##from ligand_results import ligand_results
##from receptor_results import receptor_results
##all_ligand_results = []
##all_receptor_results = []
##for stem in stems: 
##    all_ligand_results.extend(ligand_results[stem])
##    all_receptor_results.extend(receptor_results[stem])
##here results are [lig1,lig2,lig3,lig4...]
## here rows in image will correspond to ligands and columns to receptors

##import numpy.oldnumeric as Numeric
##ligand_results = Numeric.array(all_ligand_results)
##ligand_results.shape = (100,100)
##now ligand_results are #[[lig1]
                         #[lig2]
                         #[lig3]
                         #[lig4]
                         #...]
##receptor_results = Numeric.array(all_receptor_results)
##receptor_results.shape = (100,100)
##now receptor_results are #[[rec1]
                           #[rec2]
                           #[rec3]
                           #[rec4]
                           #...]
## here rows in image will correspond to receptors and columns to ligands
##mini = min(results.ravel())
###clamp upper bound to 0
##maxi = 0
## from pixelMap2D import PixelMap2D
##rw = PixelMap2D(results, mini=mini, maxi=maxi)
##rw.show_image(zero_color=(255, 255, 255))


        
