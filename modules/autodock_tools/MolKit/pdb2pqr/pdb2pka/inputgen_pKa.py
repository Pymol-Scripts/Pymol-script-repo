#!/bin/env python
#
# $Id: inputgen_pKa.py,v 1.1 2007/08/03 20:59:24 sargis Exp $
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
            defaults['coarsedim'].append(axis*2.0)
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
                #print split
                #print '0123456789012345678901234567890123456789012345678901234567890123456789'
                #print line
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
        self.type=type
        if type=='desolv':
            self.set_method('mg-manual')
        elif type=='background': 
            self.set_method('mg-auto')
        elif type=='intene':
            self.set_method('mg-auto')
            self.setfineCenter(self.coarsecent)
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
        text += "end\n"
        text += "elec\n"
        text += "    %s\n" % self.method
        text += "    dime %i %i %i\n" % (self.finegridpoints[0], self.finegridpoints[1], self.finegridpoints[2])
        if self.method=='mg-auto':
            text += "    cglen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    cgcent %.3f %.3f %.3f\n" %(self.coarsecent[0],self.coarsecent[1],self.coarsecent[2])
        
            text += "    fglen %.4f %.4f %.4f\n" % (self.finedim[0], self.finedim[1], self.finedim[2])
            text += "    fgcent %.3f %.3f %.3f\n" %(self.finecent[0],self.finecent[1],self.finecent[2])
        elif self.method=='mg-manual':
            text += "    glen %.4f %.4f %.4f\n" % (self.coarsedim[0], self.coarsedim[1], self.coarsedim[2])
            text += "    gcent %.3f %.3f %.3f\n" %(self.coarsecent[0],self.coarsecent[1],self.coarsecent[2])
            text += "    nlev 4\n"
            
        else:
            raise 'unknown method'
        #
        text += "    mol 1\n"                            
        text += "    lpbe\n"                             
        text += "    bcfl sdh\n"                           
        text += "    ion 1 0.150 2.0\n"            
        text += "    ion -1 0.150 2.0\n"           
        text += "    pdie %5.2f\n"  %self.pdie                
        text += "    sdie %5.2f\n" %self.sdie                
        text += "    srfm smol\n"                   
        text += "    chgm spl2\n"
        text += "    srad 1.4\n"          
        text += "    swin 0.3\n"         
        text += "    temp 298.15\n"     
        text += "    gamma 0.105\n"    
        text += "    calcenergy total\n"
        text += "    calcforce no\n"
        text += "    write pot dx pot\n"
        text += "    write smol dx acc\n"
        text += "end\n"
        return text

    #
    # ------
    #

    def getText(self):
        #
        # Energy statements
        #
        text=self.getText_sub()
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
        raise 'type not set'

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
