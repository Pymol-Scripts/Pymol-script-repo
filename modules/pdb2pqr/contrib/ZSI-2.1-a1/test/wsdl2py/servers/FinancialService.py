#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import sys
from FinancialService_server import *
from ZSI.ServiceContainer import AsServer

class Service(FinancialService):

    def soap_getPV(self, ps):
        request,response = FinancialService.soap_getPV(self, ps)
        args = request

        # Worker code: Actual present value calculation
        t = 0
        PV = 0.0
        for CF in args._CFSequence._CF:
            PV += (CF or 0.0) * ((args._irate / 100.0 + 1) ** (-t)) 
            t += 1

        #print "Present value is: ", PV

        # assign return values to response object
        #class SimpleTypeWrapper(float): typecode = getPVResponseWrapper()
        # WARNING specify value eg. SimpleTypeWrapper(1)
        #response = SimpleTypeWrapper(PV)

        response = response.__class__(PV)
        return request,response



if __name__ == "__main__" :
    port = int(sys.argv[1])
    AsServer(port, (Service('test'),))
