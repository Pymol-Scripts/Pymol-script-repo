#!/usr/bin/env python

############################################################################
# David W. Robertson, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys, unittest
from ServiceTest import ServiceTestCase, ServiceTestSuite

import re
from ZSI import EvaluateException
"""
Unittest for contacting the TerraService Web service.

WSDL:  http://terraservice.net/TerraService.asmx?WSDL
"""

CONFIG_FILE = 'config.txt'
CONFIG_SECTION = 'complex_types'
SERVICE_NAME = 'TerraService'
PORT_NAME = 'TerraServiceSoap'
EXCEPTION_STRING_SERIALIZE = r"Serializing ConvertPlaceToLonLatPt xmlns=\"http://terraserver-usa.com/terraserver/\"._place, Exception Serializing place xmlns=\"http://terraserver-usa.com/terraserver/\"._City, AttributeError 'int' object has no attribute \'replace\'"

SERIALIZE_PATTERN = re.compile(EXCEPTION_STRING_SERIALIZE)

class TerraServiceTest(ServiceTestCase):
    """Test case for TerraService Web service
    """
    name = "test_TerraService"

    def test_ConvertPlaceToLonLatPt(self):
        operationName = 'ConvertPlaceToLonLatPt'
        request = self.getInputMessageInstance(operationName)
        request._place = self._moduleDict[self._typeModuleName].ns1.Place_Def()
        request._place._City = 'Oak Harbor'
        request._place._State = 'Washington'
        request._place._Country = 'United States'
        response = self.RPC(operationName, request)
        


    def test_ConvertLonLatPtToNearestPlace(self):
        operationName = 'ConvertLonLatPtToNearestPlace'
        request = self.getInputMessageInstance(operationName)
        request._place = self._moduleDict[self._typeModuleName].ns1.Place_Def()
        request._point = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._point._Lon = -122.643
        request._point._Lat = 48.297
        response = self.RPC(operationName, request)   
    
    def test_ConvertLonLatPtToUtmPt(self):
        operationName = 'ConvertLonLatPtToUtmPt'
        request = self.getInputMessageInstance(operationName)
        request._point = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._point._Lon = -122.643
        request._point._Lat = 48.297
        response = self.RPC(operationName, request)  

    def test_ConvertUtmPtToLonLatPt(self):
        operationName = 'ConvertUtmPtToLonLatPt'
        request = self.getInputMessageInstance(operationName) 
        request._utm = self._moduleDict[self._typeModuleName].ns1.UtmPt_Def()
        request._utm._X =  526703.512403
        request._utm._Y =  5348595.96493
        request._utm._Zone =  10
        response = self.RPC(operationName, request)  

    def test_CountPlacesInRect(self):
        operationName = 'CountPlacesInRect'
        request = self.getInputMessageInstance(operationName)
        request._upperleft = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._upperleft._Lon = -122.647
        request._upperleft._Lat = 48.293
        request._lowerright = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._lowerright._Lon = request._upperleft._Lon + 1.0
        request._lowerright._Lat = request._upperleft._Lon - 1.0
        request._ptype = "HillMountain"
        response = self.RPC(operationName, request)
    
    def test_GetAreaFromPt(self):
        operationName = 'GetAreaFromPt'
        request = self.getInputMessageInstance(operationName)
        request._center = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._center._Lon = -122.647
        request._center._Lat = 48.293
        request._theme = 'Topo'
        request._scale = "Scale2m"
        request._displayPixWidth = 2
        request._displayPixHeight = 2
        response = self.RPC(operationName, request)

    def test_GetAreaFromRect(self):
        operationName = 'GetAreaFromRect'
        request = self.getInputMessageInstance(operationName)
        request._upperLeft = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._upperLeft._Lon = -122.647
        request._upperLeft._Lat = 48.293
        request._lowerRight = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._lowerRight._Lon = request._upperLeft._Lon + 1.0
        request._lowerRight._Lat = request._upperLeft._Lat - 1.0
        request._theme = 'Topo'
        request._scale = "Scale2m"
        response = self.RPC(operationName, request)

    def test_GetAreaFromTileId(self):
        operationName = 'GetAreaFromTileId'
        request = self.getInputMessageInstance(operationName)
        id = self._moduleDict[self._typeModuleName].ns1.TileId_Def()
        id._Theme = 'Topo'
        id._Scale = "Scale2m"
        id._Scene = 8
        id._X = 20
        id._y = 20
        request._id = id
        request._displayPixWidth = 2
        request._displayPixHeight = 2
        response = self.RPC(operationName, request)

    def test_GetLatLonMetrics(self):
        operationName = 'GetLatLonMetrics'
        request = self.getInputMessageInstance(operationName)
        request._point = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._point._Lon = -122.647
        request._point._Lat = 48.293
        response = self.RPC(operationName, request)

        # derived type (enum) problem
        # skipping it for now

              
        # derived type (enum) problem
        # also inconsistent timeout problem for this call


    def test_GetPlaceListInRect(self):
        operationName = 'GetPlaceListInRect'
        request = self.getInputMessageInstance(operationName)
        request._upperleft = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._upperleft._Lon = -123.0
        request._upperleft._Lat = 44.0
        request._lowerright = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
            # needs to be small, otherwise different items
            # returned each time
        request._lowerright._Lon = -122.8
        request._lowerright._Lat = 43.8
        request._ptype = "HillMountain"
        request._MaxItems = 3
        response = self.RPC(operationName, request)

    def test_GetTheme(self):
        operationName = 'GetTheme'
        request = self.getInputMessageInstance(operationName)
        request._theme = 'Topo'
        response = self.RPC(operationName, request)

    def test_GetTile(self):
        operationName = 'GetTile'
        request = self.getInputMessageInstance(operationName)
        request._id = self._moduleDict[self._typeModuleName].ns1.TileId_Def()
        request._id._Theme = 'Topo'
        request._id._Scale = 'Scale2m'
        request._id._Scene = 8
        request._id._X = 20
        request._id._Y = 20
        response = self.RPC(operationName, request)

    def test_GetTileMetaFromLonLatPt(self):
        operationName = 'GetTileMetaFromLonLatPt'
        request = self.getInputMessageInstance(operationName)
        request._theme = 'Topo'
        request._point = self._moduleDict[self._typeModuleName].ns1.LonLatPt_Def()
        request._point._Lon = -122.64
        request._point._Lat = 48.29
        request._scale = "Scale4m"
        response = self.RPC(operationName, request)

    def test_GetTileMetaFromTileId(self):
        operationName = 'GetTileMetaFromTileId'
        request = self.getInputMessageInstance(operationName)
        request._id = self._moduleDict[self._typeModuleName].ns1.TileId_Def()
        request._id._Theme = 'Topo'
        request._id._Scale = 'Scale2m'
        request._id._Scene = 8
        request._id._X = 20
        request._id._Y = 20
        response = self.RPC(operationName, request)


class TerraServiceTestFailures(ServiceTestCase):
    name = "test_TerraService"
    
    def test_ConvertPlaceToLonLatPt_x1(self):
        """
        This test should fail
        """
        operationName = 'ConvertPlaceToLonLatPt'
        request = self.getInputMessageInstance(operationName)
        request._place = self._moduleDict[self._typeModuleName].ns1.Place_Def()
        request._place._City = 1
        request._place._State = 'Washington'
        request._place._Country = 'United States'
        
        try:
            response = self.RPC(operationName, request)
            
        except Exception, msg:
            exceptionString = str(msg)
            if SERIALIZE_PATTERN.match(exceptionString):
                pass
            else:
                raise

    def test_GetPlaceFacts(self):
        operationName = 'GetPlaceFacts'
        request = self.getInputMessageInstance(operationName)
        request._place = self._moduleDict[self._typeModuleName].ns1.Place_Def()
        request._place._City = 'Seattle'
        request._place._State = 'Washington'
        request._place._Country = 'United States'
        try:
            response = self.RPC(operationName, request)
        except EvaluateException, ex:
            pass

    def test_GetPlaceList(self):
        operationName = 'GetPlaceList'
        request = self.getInputMessageInstance(operationName)
        request._placeName = 'New York'
        request._MaxItems = 5
        request._imagePresence = 0
        try:
            response = self.RPC(operationName, request)
        except EvaluateException, ex:
            pass

def makeTestSuite():
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TerraServiceTest, 'test_'))
    suite.addTest(unittest.makeSuite(TerraServiceTestFailures, 'test_'))
    return suite


if __name__ == "__main__" :
    unittest.TestProgram(defaultTest="makeTestSuite")
