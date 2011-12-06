#
#
#

from EchoServer_client import *
import sys, time

TRACE=None
loc = EchoServerLocator()
port = loc.getEchoServer(url='http://localhost:8000/echo', tracefile=TRACE)

msg = EchoRequest()
msg.Value = 1
rsp = port.Echo(msg)
print "INTEGER: ", rsp.Value

msg.Value = "HI"
rsp = port.Echo(msg)
print "STRING: ", rsp.Value

msg.Value = 1.10000
rsp = port.Echo(msg)
print "FLOAT: ", rsp.Value


msg.Value = dict(milk=dict(cost=3.15, unit="gal"))
rsp = port.Echo(msg)
print "DICT: ", rsp.Value
