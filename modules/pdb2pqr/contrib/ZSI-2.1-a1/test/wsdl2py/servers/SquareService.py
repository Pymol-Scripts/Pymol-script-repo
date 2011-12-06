#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys
import SquareService_server as Square
from ZSI.ServiceContainer import AsServer

class Service(Square.SquareService):

    def soap_getSquare(self, ps):
        request,response = Square.SquareService.soap_getSquare(self, ps)
        response._return = self.getSquare(request._x)
        return request,response

    def getSquare(self, x):
        return x**2


if __name__ == "__main__" :
    port = int(sys.argv[1])
    AsServer(port, (Service('test'),))
