## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
#
# Author: Sophie COON, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################


# $Header: /opt/cvs/python/packages/share1.5/mglutil/util/colorUtil.py,v 1.15 2007/07/24 17:30:40 vareille Exp $
#
# $Id: colorUtil.py,v 1.15 2007/07/24 17:30:40 vareille Exp $
#

import numpy.oldnumeric as Numeric, math
from types import FunctionType

""" This Python module implements a set of methods and classes to manipulate
colors representation.
Should put ColorMap here later on...
"""

def TkColor(col):
    """
    col can be a :
    - RGB triplet of int (0-255) or floats (0.0, 1.0)
    """
    if max(col)<=1.0: col = map( lambda x: x*255, col)
    return '#%02X%02X%02X' % (col[0],col[1],col[2])

def ToHEX(col, mode="RGB", flag255=0):
    """
    ToHEX takes a :
    - RGB triplet of int (0-255) or floats (0.0, 1.0)
    - HSV triplet of float (0.0, 1.0)
    and converts it to a string representing the hexadecimal value
    '#%02X%02X%02X'
    """
    # If HSV need to convert it to RGB
    assert mode in ['RGB','HSV']
    assert flag255 in [0,1]
    if mode == "HSV":
        col = ToRGB(col)
    # If float (0.0-1.0) need to transform it into (0-255)
    if flag255 == 0:
        col = Numeric.array(col,'f')*255
        col.tolist()
    return '#%02X%02X%02X'%(col[0],col[1],col[2])

def ToRGBA(col, mode="HSV", flag255=0):
    assert mode in ['HSV']
    if mode=="HSV":
        assert len(col)==4
        rgbCol = list(col)
        rgb = ToRGB(col[:3],mode, flag255)
        if not rgb is None:
            rgbCol[:3] = ToRGB(col[:3],mode, flag255)
        return rgbCol
    elif mode == 'HEX':
        print 'I am not sure of what to do ?'
    
def ToRGB(col, mode="HSV", flag255=0):
    """
    ToRGB takes a:
    - HSV triplet
    - HEX string
    and returns a the corresponding rgb triplet (0.0 to 1.0) or (0-255) if
    the 255Flag is set to 1.
    """
    assert mode in ['HSV', 'HEX']
    if mode == 'HEX':
        if not col[0] == '#' or len(col)!=7:
            print 'invalid HEX color needs to be "#rrggbb" '
        else:
            r = float(eval('0x'+col[1:3]))
            g = float(eval('0x'+col[3:5]))
            b = float(eval('0x'+col[5:]))
            #print a
            if flag255 == 0:
                return (r/255.,g/255.,b/255.)
            else:
                return (int(r), int(g), int(b))
            
    elif mode == 'HSV':
        l = len(col)
        assert l==3
        assert max(col) <= 1.0
        assert min(col) >= 0.0
        
        v = col[2]
        if v == 0.0:
            return (0.0, 0.0, 0.0)

        s = col[1]
        if s == 0.0:
            if flag255:
                nCol = Numeric.array((v,v,v),'f')*255
                return tuple(nCol)
            return (v, v, v)

        h = col[0]*6.0
        if h>=6.0: h = 0.0
        i = int(h)
        f = h - i
        p = v*(1.0 - s)
        q = v*(1.0-(s*f))
        t = v*(1.0-s*(1.0-f))

        if i==0:
            if flag255:
                nCol = Numeric.array((v,t,p),'f')*255
                return tuple(nCol)
            return (v,t,p)
        elif i==1:
            if flag255:
                nCol = Numeric.array((q,v,p),'f')*255
                return tuple(nCol)
            return (q,v,p)
        elif i==2:
            if flag255:
                nCol = Numeric.array((p,v,t),'f')*255
                return tuple(nCol)
            return (p,v,t)
        elif i==3:
            if flag255:
                nCol = Numeric.array((p,q,v),'f')*255
                return tuple(nCol)
            return (p,q,v)
        elif i==4:
            if flag255:
                nCol = Numeric.array((t,p,v),'f')*255
                return tuple(nCol)
            return (t,p,v)
        elif i==5:
            if flag255:
                nCol = Numeric.array((v,p,q),'f')*255
                return tuple(nCol)
            return (v,p,q)
        else:
            print "botch in col_to_rgb"


def Hue2RGB( v1, v2, vH ):
    """ used by fromHslToRgb
"""
    if vH < 0:
        vH += 1
    if vH > 1:
        vH -= 1
    if ( 6 * vH ) < 1:
        return v1 + ( v2 - v1 ) * 6 * vH
    if ( 2 * vH ) < 1:
        return v2
    if ( 3 * vH ) < 2:
        return v1 + ( v2 - v1 ) * ( ( 2 / 3. ) - vH ) * 6
    return v1


def HSL2RGB(H,S,L):
    """
HSL values = 0 - 1
RGB results = 0 - 1
"""
    if S == 0: 
       R = L
       G = L
       B = L
    else:
       if L < 0.5 :
           var_2 = L * ( 1 + S )
       else:
           var_2 = ( L + S ) - ( S * L )
       var_1 = 2 * L - var_2
       R = Hue2RGB( var_1, var_2, H + ( 1 / 3. ) )
       G = Hue2RGB( var_1, var_2, H )
       B = Hue2RGB( var_1, var_2, H - ( 1 / 3. ) )
    return [R,G,B]


def RGB2HSL(R,G,B):
    """
RGB values = 0 - 1
HSL results = 0 - 1
"""   
    var_Min = min( R, G, B )    
    var_Max = max( R, G, B )    
    del_Max = var_Max - var_Min             
    L = ( var_Max + var_Min ) / 2
    if del_Max == 0: 
       H = 0                               
       S = 0
    else:
       if L < .5:
           S = del_Max / ( var_Max + var_Min )
       else:
           S = del_Max / ( 2. - var_Max - var_Min )
       del_R = ( ( ( var_Max - R ) / 6. ) + ( del_Max / 2. ) ) / del_Max
       del_G = ( ( ( var_Max - G ) / 6. ) + ( del_Max / 2. ) ) / del_Max
       del_B = ( ( ( var_Max - B ) / 6. ) + ( del_Max / 2. ) ) / del_Max
       if R == var_Max:
           H = del_B - del_G
       elif G == var_Max:
           H = ( 1 / 3. ) + del_R - del_B
       elif B == var_Max:
           H = ( 2 / 3. ) + del_G - del_R
       if H < 0:
           H += 1
       if H > 1:
            H -= 1
    return [H, S, L]


def RGB2HSL_list(RGB):
    """
RGB values = 0 - 1
HSL results = 0 - 1
"""
    return RGB2HSL(RGB[0], RGB[1], RGB[2])


def RGBA2HSLA_list(RGBA):
    """
RGB values = 0 - 1
HSL results = 0 - 1
"""
    hsla = RGB2HSL(RGBA[0], RGBA[1], RGBA[2])
    hsla.append(RGBA[3])
    return hsla


def HSLA2RGBA_list(HSLA):
    """
RGB values = 0 - 1
HSL results = 0 - 1
"""
    rgba = HSL2RGB(HSLA[0], HSLA[1], HSLA[2])
    rgba.append(HSLA[3])
    return rgba


def ToHSV(col, mode="RGB", flag255=0):
    """
ToHSV takes a:
- RGB triplet (0-255 or 0.0-1.0)
- HEX string ('#%02x%02x%02x'%(rr,gg,bb))
and returns a the corresponding HSV triplet (0.0 to 1.0)
"""
    assert mode in ['HEX','RGB']
    assert flag255 in [0,1]
    if mode == "HEX":
        if col[0] != '#' or len(col)!=7:
            print 'invalid HEX color needs to be "#rrggbb" '
        else:
            r = float(eval('0x'+col[1:3]))
            g = float(eval('0x'+col[3:5]))
            b = float(eval('0x'+col[5:]))
            newCol = Numeric.array([r,g,b], 'f')
            assert min(newCol) >=0.0
            assert max(newCol) <=255.0
            col= Numeric.array(newCol,'f')/255
    elif mode == "RGB":
        assert len(col)==3
        if flag255: col =  Numeric.array(col,'f')/255
    else:
        print "color mode should be 'HEX' or 'RGB'"
    maxi = max(col)
    mini = min(col)
    r,g,b = col
    assert maxi<= 1.0 
    assert mini >= 0.0
    if maxi > 0.0001: s = (maxi - mini)/maxi
    else: s = 0.0
    if s < 0.0001: h = 0.0
    else:
        delta = maxi - mini
        if r == maxi: h = (g - b)/delta
        elif g == maxi: h = 2.0 + (b - r)/delta
        elif b == maxi: h = 4.0 + (r - g)/delta
        h = h/6.0
        if h < 0.0: h = h + 1.0

    return (h,s,maxi)


def ToHSVA(rgb):
    pass


#class ColorPalette0:
#    FLOAT = 0
#    INT = 1
#
#    def __init__(self, name, colorDict={}, readonly=0, colortype=None,
#                 info='', sortedkeys=None, lookupMember=None):
#        self.name = name
#        self.readonly = readonly
#        self.colors = colorDict
#        self.info = info
#        self.viewer = None
#        self.sortedkeys = sortedkeys
#        if colortype is None:
#            self.colortype = self.FLOAT
#        self.lookupMember = lookupMember
#
#    def _lookup(self, name):
#        if not name in self.colors.keys():
#            return (0., 1., 0.)
#        return self.colors[name]
#
#    def lookup(self, objects):
#        # Maybe should try that first in case all the objects don't have the
#        # lookup member
#        names = objects.getAll(self.lookupMember)
#        return map( self._lookup, names)
#
#    def display(self,*args, **kw):
#        """ Will create an instance of PaletteChooser later on"""
#        pass
#    
#    def undisplay(self, *args, **kw):
#        pass
#
#    def copy(self):
#        """make a deep copy of a palette"""
#        import copy
#        c = copy.copy(self)
#        c.readonly = 0
#        c.colors = copy.deepcopy(self.colors)
#        if self.sortedkeys:
#            c.sortedkeys = self.sortedkeys[:]
#        return c
#
#
#class ColorPaletteFunction0(ColorPalette0):
#    def __init__(self, name, colorDict={}, readonly=0, colortype=None,
#                 info='', sortedkeys=None, lookupFunction = None):
#        """ lookupFunction : needs to be function or a lambda function"""
#        ColorPalette.__init__(self, name, colorDict, readonly,colortype,
#                               info, sortedkeys)
#        if not type(lookupFunction) is FunctionType:
#            self.lookupFunction = None
#
#        self.lookupFunction = lookupFunction
#                               
#    def lookup(self, objects):
#        # maybe should do that in a try to catch the exception in case it
#        # doesnt work
#        names = map(self.lookupFunction, objects)
#        return map(self._lookup, names)
