#
#
#
#$Id: test_colorUtil.py,v 1.2 2005/06/17 16:23:49 sowjanya Exp $

##############################################################################
##  
##  Authors : Sowjanya Karnati, Michel F Sanner
##
##############################################################################



import sys,unittest,Tkinter
from time import sleep
from mglutil.util.colorUtil import *

class ColorUtilBaseTest(unittest.TestCase):

#####TkColor

    def test_TkColor_1(self):
        """tests TkColor with RGB float triplet"""
        col=(0.0,1.0,0.0)
        self.rvalue=TkColor(col)
        self.assertEqual(self.rvalue,'#00FF00')
        

    def test_TkColor_2(self):
        """tests TkColor with inavalid col arguement"""
        col=(0.0,1.0)
        self.assertRaises(IndexError,TkColor,col)


    def test_TkColr_3(self):
        """tests TkColor with inavalid col arguement"""
        col ='hello'
        self.assertRaises(TypeError,TkColor,col)

    def test_TkColor_4(self):
        """tests TkColor with RGB Triplet"""
        col = (240,248,255)
        self.rvalue=TkColor(col)
        self.assertEqual(self.rvalue,"#F0F8FF")

######ToHEX

    def test_ToHEX_1(self):
        """tests ToHEX with mode=HSV"""
        color = (0.0,1.0,0.5)
        self.rvalue =ToHEX(color,mode='HSV')
        self.assertEqual(self.rvalue,'#7F0000')
        
    def test_ToHEX_2(self):
        """tests ToHEX with mode=RGB"""
        color=(0.0,1.0,0.5)
        self.rvalue =ToHEX(color,mode='RGB')
        self.assertEqual(self.rvalue,'#00FF7F')

    def test_ToHEX_3(self):
        """tests ToHEX with flag,mode RGB"""
        color=(0.0,1.0,0.5)
        self.rvalue =ToHEX(color,mode='RGB',flag255=1)
        self.assertEqual(self.rvalue,'#000100')
    
    def test_ToHEX_4(self):
        """tests ToHEX with flag,mode HSV"""
        color=(0.0,1.0,1.0)
        self.rvalue =ToHEX(color,mode='HSV',flag255=1)
        self.assertEqual(self.rvalue,'#010000')

    def test_ToHEX_5(self):
        """tests ToHEX with invalid color arguement"""
        color=(0.0,1.0)
        self.assertRaises(IndexError,ToHEX,color)

    def test_ToHEX_7(self):
        """tests ToHEX with invalid mode arguement"""
        color=(0.0,1.0,0.5)
        mode="hello"
        self.assertRaises(AssertionError,ToHEX,color,mode)

    
    def test_ToHEX_8(self):
        """tests ToHex with invalid flag arguement"""
        color=(0.0,1.0,0.5)
        mode="HSV"
        flag255="hello"
        self.assertRaises(AssertionError,ToHEX,color,mode,flag255)

    def test_ToHex_9(self):
        """tests ToHex with RGB color Triplet"""
        col = (240,248,255)
        self.rvalue =ToHEX(col,mode='RGB')

######ToRGBA

    def test_ToRGBA_1(self):
        """tests ToRGBA """
        color=(0.0,1.0,0.5,0.5)
        self.rvalue =ToRGBA(color) 
        self.assertEqual(self.rvalue,[0.5, 0.0, 0.0, 0.5])
        
    def test_ToRGBA_2(self):
        """Tests ToRGBA with invalid mode"""
        color=(0.0,1.0,0.5,0.5)
        self.assertRaises(AssertionError,ToRGBA,color,mode='ABCD')        


    def test_ToRGBA_3(self):
        """tests ToRGBA with flag255"""
        color=(0.0,1.0,0.5,0.5)
        self.rvalue =ToRGBA(color,flag255=1) 
        self.assertEqual(self.rvalue,[127.5, 0.0, 0.0, 0.5])
    
   
########ToRGB

    def test_ToRGB_1(self):
        """tests ToRGB """
        color=(0.0,1.0,0.5)
        self.rvalue = ToRGB(color)
        self.assertEqual(self.rvalue,(0.5, 0.0, 0.0))

    def test_ToRGB_2(self):
        """tests ToRGB with mode HEX"""
        color='#FFFF00'
        self.rvalue = ToRGB(color,mode='HEX')    
        self.assertEqual(self.rvalue,(1.0, 1.0, 0.0))

    def test_ToRGB_3(self):
        """tests ToRGB with invalid color"""
        color="hello"
        self.assertRaises(AssertionError,ToRGB,color)
        
    def test_ToRGB_4(self):
        """tests ToRGB with flag"""
        color=(0.0,1.0,0.5)
        self.rvalue = ToRGB(color,flag255=1)
        self.assertEqual(self.rvalue,(127.5, 0.0, 0.0))

    def test_ToRGB_5(self):
        """tests ToRGB with mode HEX anf flag255 1"""
        color='#FF0000'
        self.rvalue = ToRGB(color,mode='HEX',flag255=1)
        self.assertEqual(self.rvalue,(255, 0, 0))

    def test_ToRGB_6(self):
        """tests ToRGB with invalid mode"""
        color='#FF0000'
        self.assertRaises(AssertionError,ToRGB,color,mode='hello')



#####ToHSV

    def test_ToHSV_1(self):
        """tests ToHSV """
        color=(0.0,1.0,0.5)
        self.rvalue = ToHSV(color)
        self.assertEqual(self.rvalue,(0.41666666666666669, 1.0, 1.0))

    def test_ToHSV_2(self):
        """tests ToHSV with mode HEX"""
        self.rvalue = ToHSV(col='#FFFF00',mode='HEX')    
        self.assertEqual(self.rvalue,((0.16666666666666666, 1.0, 1.0)))

    def test_ToHSV_3(self):
        """tests ToHSV with invalid color"""
        color="hello"
        self.assertRaises(AssertionError,ToHSV,color)
      
    
    def test_ToHSV_5(self):
        """tests ToHSV with flag255 in valid """
        color='#FF0000'
        self.assertRaises(AssertionError,ToHSV,color,mode="HEX",flag255="hello")

    def test_ToHSV_6(self):
        """tests ToHSV with mode invalid"""
        color='#FF0000'
        self.assertRaises(AssertionError,ToHSV,color,mode='hello')       

    def test_ToHSV_7(self):
        """tests ToHSV with RGB color"""
        col = (0,255,0)
        self.rvalue = ToHSV(col,flag255=1)
        self.assertEqual(self.rvalue,(0.33333333333333331, 1, 1))

    
 
if __name__ == '__main__':
   unittest.main() 


        
