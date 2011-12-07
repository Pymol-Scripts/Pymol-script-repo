## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#
# Author Daniel Stoffler  (Dec 2001)          Copyright D. Stoffler, TSRI

import numpy.oldnumeric as Numeric, string

class VRML2IndexedFaceSetToCoords:
    """ VRML2 IndexedFaceSet -> coords & indices.
    Please note this is BY NO MEANS a VRML parser but rather a hack that
    allows to convert SIMPLE VRML2 indexedFaceSets into coords and indices
    which can be displayed using the IndexedPolygon node of the  ViPEr
    environment """
    
    def init(self, data=None):
        self.data = data
              

    def set(self, data):
        self.data = data


    def compute(self):
        spl=[]
        coordflag = 0
        indexflag = 0
        scaleflag = 0
        coords = []
        indices = []
        scale=[]
        llc = []
        lli = []
        result = []
        ind = []
        length = 0
        
        for i in range(len(self.data)):
            d = self.data[i]
            spl = string.split(d)

            for k in range(len(spl)):
                if spl[k] == 'coord' and spl[k+1] == 'Coordinate':
                    coordflag = 1

                if spl[k] == 'coordIndex':
                    indexflag = 1

                if spl[k] == 'scale':
                    scaleflag = 1



            if coordflag == 1:
                if spl[0]=='}': # we reach the end of 'coords'
                    coordflag=0
                    continue

                if spl[0] == 'coord':  #this is the first line which we dont
                    continue           #want to parse

                if spl[-1] == ']':     #this is the last line where we want
                    lc = spl[-4:-1]    #to get rid of the ']'
                else:
                    lc = spl[-3:]
                    lc[-1] = lc[-1][:-1]

                llc=[]
                for n in range(len(lc)):

                    llc.append(float(lc[n]))
                coords.append(llc)


            if indexflag == 1:
                testEnd = string.split(d)
#                if testEnd[0]=='texCoord':
#                    indexflag=0
#                    continue

                if spl[-1] == ']':      #this is the last line where we want
                    li = spl[-9:-1]     #to get rid of the ']'
                    li[-1] = li[-1]+',' #and add a ',' to the last number
                    indexflag=0
                else:
                    li = spl[-8:]

                ind.extend(li)
                lli = []


            if scaleflag == 1:
                sc = string.split(d)
                scale.append( float(sc[1]) )
                scale.append( float(sc[2]) )
                scale.append( float(sc[3]) )
                scaleflag = 0

        # make index list. Note that we add -1 as many times as needed
        # to each index entry in order to make them all the same length
        
        lenght=0 # this is the variable that tells us how many -1 we will ad
        for n in range(len(ind)):
            if ind[n] != '-1,':
                lli.append(int(ind[n][:-1]))
            else:
                lli.reverse() # index list has to be inverted
                lli.append(-1)
                if len(lli) > length: length = len(lli)
                indices.append(lli)
                lli = []



        # here we go again over the indices and add the -1 if necessary
        for m in range(len(indices)):
            if len(indices[m]) < length:
                for o in range(length-len(indices[m])):
                    indices[m].append(-1)
                
        # apply scale to coords
        if len(scale):
            coords = Numeric.multiply(coords, scale)
        
        result = [coords]+[indices]
        return result
    

    def __call__(self):
        out=[]
        out = self.compute()
        return out
