#!/usr/bin/env python
############################################################################
# Joshua R. Boverhof, LBNL
# See LBNLCopyright for copyright notice!
###########################################################################
import os, sys, unittest
from ServiceTest import main, ServiceTestCase, ServiceTestSuite
from ZSI import FaultException, ParsedSoap, SoapWriter
"""
Unittest 

WSDL:   
"""

# General targets
def dispatch():
    """Run all dispatch tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_dispatch'))
    return suite

def local():
    """Run all local tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_local'))
    return suite

def net():
    """Run all network tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_net'))
    return suite
    
def all():
    """Run all tests"""
    suite = ServiceTestSuite()
    suite.addTest(unittest.makeSuite(TestCase, 'test_'))
    return suite


class TestCase(ServiceTestCase):
    name = "test_Racing"
    client_file_name = "Racing_client.py"
    types_file_name  = "Racing_types.py"
    server_file_name = "Racing_server.py"

    def __init__(self, methodName):
        ServiceTestCase.__init__(self, methodName)
        self.wsdl2py_args.append('-b')

    def test_local_anyType(self):
        """rpc/lit, testing if <any/> lax content handling
        should get back dicts and strings 
        """
        ps = ParsedSoap(MSG)
        pyobj = ps.Parse(self.client_module.EventApproximatesSoapOut.typecode)

        any = {'PoolTotals': {'Pool': {'Total': u'4117.66', 'ENumbers': None, 'JackpotNet': None}}, 'Approximates': {'Pool': {'Win': u'3.90,0.00,10.40,11.80,4.70,29.50,29.90,2.40,19.80,0.00', 'Place': u'1.04,0.00,2.80,5.90,2.00,5.20,7.40,1.04,4.00,0.00'}}}

        self.failUnless(pyobj.EventApproximatesResult.Any == any, 'Failed match:\n %s\n\n%s' %(
            pyobj.EventApproximatesResult.Any, any))


        pyobj.EventApproximatesResult.Any = dict(pyobj.EventApproximatesResult.Any)
        sw = SoapWriter()
        sw.serialize(pyobj)
        print str(sw)
        ps2 = ParsedSoap(str(sw))
        pyobj2 = ps.Parse(self.client_module.EventApproximatesSoapOut.typecode)
        print "EAR: ", pyobj2.EventApproximatesResult
        print "Any: ", pyobj2.EventApproximatesResult.Any
        
        
        self.failUnless(pyobj.EventApproximatesResult.Any == pyobj2.EventApproximatesResult.Any,
            'Failed match:\n %s\n\n%s' %(pyobj.EventApproximatesResult.Any, pyobj2.EventApproximatesResult.Any))



MSG="""<?xml version="1.0" encoding="utf-8"?>
<soap:Envelope xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema">
	<soap:Body>
		<EventApproximatesResponse xmlns="http://direct.tab.com.au/LiveOdds/">
			<EventApproximatesResult>
				<Information xmlns="" Jurisdiction="NSW">
					<Approximates RID="SG_20070123_05" Timestamp="20:17:17 20070123">
						<Pool>
							<Win>3.90,0.00,10.40,11.80,4.70,29.50,29.90,2.40,19.80,0.00</Win>
							<Place status="">1.04,0.00,2.80,5.90,2.00,5.20,7.40,1.04,4.00,0.00</Place>
						</Pool>
					</Approximates>
					<PoolTotals RID="SG_20070123_05">
						<Pool BetTypeDesc="Exacta">
							<Total>451.00</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="First Four">
							<Total>1001.00</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="Place">
							<Total>1750.59</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="Quaddie">
							<Total>6052.50</Total>
							<ENumbers>05 06 07 08</ENumbers>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="Quinella">
							<Total>865.00</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="Running Double">
							<Total>21.50</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="Trifecta">
							<Total>2575.50</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
						<Pool BetTypeDesc="Win">
							<Total>4117.66</Total>
							<ENumbers/>
							<JackpotNet/>
						</Pool>
					</PoolTotals>
				</Information>
			</EventApproximatesResult>
		</EventApproximatesResponse>
	</soap:Body>
</soap:Envelope>"""


if __name__ == "__main__" :
    main()
