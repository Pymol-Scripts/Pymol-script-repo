#!/bin/env python
#
# $Id: inputgen_pKa.py 974 2011-06-30 14:58:03Z jens_nielsen $
#


class inputGen:
    #
    # Creating input files for the pKa calculations
    #
    
    def __init__(self, pqrpath):
        """
            Initialize the inputGen class
            
            Parameters
                pqrpath: The full path to the PQR file (string)
                center_residue: The residue on which the fine grid will be centered
        """
        #
        # Defs by Jens
        #
        defaults={'sdie':78.54,'pdie':8.0,
                  'finegridpoints':[65,65,65],
                  'coarsedim':[],
                  'method':'mg-auto'}

        #
        # PQR file
        #
        self.pqrfile=pqrpath
        self.pqrname=pqrpath
        self.fullpath=pqrpath
        self.type='not set'
        #
        center,extent=self.getCenter()
        #
        # Make the coarse grid twice as big as the protein
        #
        for axis in extent:
            defaults['coarsedim'].append(axis*3.0)
        #
        # Center coarse grid on the center of the molecule
        #
        defaults['coarsecent']=center
        #
        # fine grid is 16.25 (res = 0.25 A/grid point)
        #
        defaults['finedim']=[16.25,16.25,16.25]
        #
        # Center the fine grid on the center of the protein to start with
        #
        defaults['finecent']=center
        #
        # Set all the attributes
        #
        for prop in defaults.keys():
            setattr(self,prop,defaults[prop])
        return

    #
    # -----
    #
    
    def getCenter(self):
        #
        # Reads the PQR file and get the dimensions and the center of the molecule
        #
        import string
        coords=[]
        fd=open(self.pqrfile)
        line=fd.readline()
        while line:
            split=string.split(line)
            if split[0] in ['ATOM','HETATM']:
                x=float(line[30:38])
                y=float(line[39:46])
                z=float(line[47:55])
                coords.append([x,y,z])
            line=fd.readline()
        fd.close()
        #
        # Find extent 
        #
        minmax=[[999.9,-999.9],[999.9,-999.9],[999.9,-999.9]]
        for coord in coords:
            for axis in xrange(3):
                if coord[axis]<minmax[axis][0]:
                    minmax[axis][0]=coord[axis]
                if coord[axis]>minmax[axis][1]:
                    minmax[axis][1]=coord[axis]
        #
        # Calc the geometric center and extent
        #
        center=[0,0,0]
        extent=[0,0,0]
        for axis in xrange(3):
            extent[axis]=minmax[axis][1]-minmax[axis][0]
            center[axis]=extent[axis]/2.0+minmax[axis][0]
        return center,extent
    
    def setfineCenter(self, center):
        #
        # Set the center of the fine grid to center
        #
        self.finecent = center
        return

    def set_method(self,method):
        #
        # Set the method
        #
        self.method=method
        return

    def set_type(self,type):
        #
        # Set the type of calculation
        #
#        print 'type: %s' % (type)
        self.type=type
        if type=='desolv':
            self.set_method('mg-manual')
        elif type=='background': 
            self.set_method('mg-auto')
            self.finedim=[self.coarsedim[0]/1.5,self.coarsedim[1]/1.5,self.coarsedim[2]/1.5]
        elif type=='intene':
            self.set_method('mg-auto')
            #self.setfineCenter(self.coarsecent)
            #
            # Set the grid to be a little bigger than the protein
            #
            self.finedim=[self.coarsedim[0]/1.5,self.coarsedim[1]/1.5,self.coarsedim[2]/1.5]
        else:
            #
            # Not a known type
            #
            raise 'unknown type',type
        return
        
    def getText_sub(self):
        """
            Get the text associated with the inputgen object

            Returns
                text:  The input file (string)
        """
        
        text  = "read\n"
        text += "    mol pqr %s\n" % self.pqrname
        if self.maps==1:
            text += "    diel dx %s %s %s\n" %(self.xdiel,self.ydiel,self.zdiel)
            text += "    kappa dx %s\n" % self.kappa
        text += "end\n"
        text += "elec\n"
        text += "    %s\n" % self.method
        text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
        if self.method=='mg-auto':
            text += "    cglen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
            text += "    cgcent %.3f %.3f %.3f\n" %(self.coarsecent[0],self.coarsecent[1],self.coarsecent[2])
            for i in range(3):
                while (self.finecent[i] - 0.5*self.finedim[i]) <= (self.coarsecent[i] - 0.5*self.coarsedim[i]):
                    self.finecent[i] += 0.1     # make sure finest mesh does not fall off the coarser meshes
                while (self.finecent[i] + 0.5*self.finedim[i]) >= (self.coarsecent[i] + 0.5*self.coarsedim[i]):
                    self.finecent[i] -= 0.1     # make sure finest mesh does not fall off the coarser meshes
            text += "    fgcent %.3f %.3f %.3f\n" %(self.finecent[0],self.finecent[1],self.finecent[2])
        elif self.method=='mg-manual':
            text += "    glen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    gcent %.3f %.3f %.3f\n" %(self.coarsecent[0],self.coarsecent[1],self.coarsecent[2])
            
        else:
            raise 'unknown method'
        #
        text += "    mol 1\n"                            
        text += "    lpbe\n"                             
        text += "    bcfl sdh\n"                           
        text += "    ion charge 1 conc 0.150 radius 2.0\n"            
        text += "    ion charge -1 conc 0.150 radius 2.0\n"           
        text += "    pdie %5.2f\n"  %self.pdie                
        text += "    sdie %5.2f\n" %self.sdie
        if self.maps == 2:
            text += "    srfm mol\n"
        else:
            text += "    srfm smol\n"
        text += "    chgm spl2\n"
        text += "    srad 1.4\n"
        text += "    swin 0.3\n" 
        text += "    sdens 10.0\n"
        text += "    temp 298.15\n"     
        text += "    calcenergy total\n"
        text += "    calcforce no\n"
        if self.maps==1:
            text += "    usemap diel 1\n"
            text += "    usemap kappa 1\n"
        text += "    write pot dx pot\n"
        text += "    write smol dx acc\n"
        if self.maps==2:
            text += "    write dielx dx xdiel_default\n"
            text += "    write diely dx ydiel_default\n"
            text += "    write dielz dx zdiel_default\n"
            text += "    write kappa dx kappa_default\n"
        text += "end\n"
        return text

    #
    # ------
    #

    def getText_sub_focus(self):
        """
            Get the text associated with the inputgen object

            Returns
                text:  The input file (string)
        """
        import math
        grids_per_A = 2.0
        dimension = 65
        scale = grids_per_A/((dimension-1.0)/(2.0*max(self.coarsedim)/3.0))
        if scale <= 1.0:
            depth = 1
        else:
            depth = int(math.log(scale)/math.log(2.0)) + 3
        
        grid_dim = float(math.pow(2,depth-1)/grids_per_A)
        
        text  = "read\n"
        text += "    mol pqr %s\n" % self.pqrname
        text += "end\n"
        text += "elec\n"
        text += "    mg-manual\n"
        text += "    dime %i %i %i\n" % (dimension, dimension, dimension)
        text += "    grid %.2f %.2f %.2f\n" % (grid_dim, grid_dim, grid_dim)
        text += "    gcent %.3f %.3f %.3f\n" % (self.finecent[0],self.finecent[1],self.finecent[2])
        text += "    mol 1\n"                            
        text += "    lpbe\n"                             
        text += "    bcfl sdh\n"                           
        text += "    ion charge 1 conc 0.150 radius 2.0\n"            
        text += "    ion charge -1 conc 0.150 radius 2.0\n"           
        text += "    pdie %5.2f\n"  %self.pdie                
        text += "    sdie %5.2f\n" %self.sdie
        text += "    srfm mol\n"
        text += "    chgm spl0\n"
        text += "    srad 1.4\n"
        text += "    swin 0.3\n" 
        text += "    sdens 10.0\n"
        text += "    temp 298.15\n"     
        text += "    calcenergy total\n"
        text += "    calcforce no\n"
        text += "    write pot dx potential0\n"
        text += "end\n"
        
        for i in range(1, depth):
            text += "elec\n"
            text += "    mg-manual\n"
            text += "    dime %i %i %i\n" % (dimension, dimension, dimension)
            text += "    grid %.2f %.2f %.2f\n" % (grid_dim, grid_dim, grid_dim)
            text += "    gcent %.3f %.3f %.3f\n" % (self.finecent[0],self.finecent[1],self.finecent[2])
            text += "    mol 1\n"                            
            text += "    lpbe\n"
            text += "    bcfl focus\n"                           
            text += "    ion charge 1 conc 0.150 radius 2.0\n"            
            text += "    ion charge -1 conc 0.150 radius 2.0\n"           
            text += "    pdie %5.2f\n"  %self.pdie                
            text += "    sdie %5.2f\n" %self.sdie
            text += "    srfm mol\n"
            text += "    chgm spl0\n"
            text += "    srad 1.4\n"
            text += "    swin 0.3\n" 
            text += "    sdens 10.0\n"
            text += "    temp 298.15\n"     
            text += "    calcenergy total\n"
            text += "    calcforce no\n"
            text += "    write pot dx potential%d\n" %i
            text += "end\n"
        return text

    #
    # ------
    #

    def getText(self):
        #
        # Energy statements
        #
        if (self.maps == 1) or (self.maps ==2):
            text=self.getText_sub()
            if self.type=='intene':
                self.setfineCenter(self.coarsecent)
        else:
            text=self.getText_sub_focus()
        if self.type=='background' or self.type=='intene':
            text += "\nprint energy 1 end\n"
            text += "\nquit\n"
            return text
        elif self.type=='desolv':
            #
            # For desolvation calcs, we calculate again with pdie=sdie
            #
            #oldpdie=self.pdie
            #self.pdie=self.sdie
            #text=text+self.getText_sub()
            #self.pdie=oldpdie
            #text += "\nprint energy 1 - 2 end\n"
            text += "\nquit\n"
            return text
        #
        # Eh?
        #
        raise Exception('type not set')

    #
    # ------
    #
    
    def printInput(self):
        import string
        period = string.find(self.fullpath,".")
        if period > 0:
            outname = self.fullpath[0:period] + ".in"
        else:
            outname = self.fullpath + ".in"
        file = open(outname, "w")
        file.write(self.getText())
        file.close()
        return outname
        
def main():
    import sys
    X=inputGenpKa(sys.argv[1],13)


if __name__ == "__main__":
    main()
