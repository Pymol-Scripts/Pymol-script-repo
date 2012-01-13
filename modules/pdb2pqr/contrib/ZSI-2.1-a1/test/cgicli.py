#! /usr/bin/env python
import getopt, socket, sys

try:
    (opts, args) = getopt.getopt(sys.argv[1:],
                    'h:p:s',
                    ( 'host=', 'port=',
                        'statusonly', 'help'))
except getopt.GetoptError, e:
    print >>sys.stderr, sys.argv[0] + ': ' + str(e)
    sys.exit(1)
if args:
    print sys.argv[0] + ': Usage error; try --help.'
    sys.exit(1)

hostname, portnum, verbose = 'localhost', 80, 1
for opt, val in opts:
    if opt in [ '--help' ]:
        print '''Options include:
    --host HOST (-h HOST)       Name of server host
    --port PORT (-p PORT)       Port server is listening on
    --statusonly (-s)           Do not output reply packets; just status code
Default is -h%s -p%d -t%s''' % \
    (hostname, portnum, ','.join([str(x) for x in tests]))
        sys.exit(0)
    if opt in [ '-h', '--host' ]:
        hostname = val
    elif opt in [ '-p', '--port' ]:
        portnum = int(val)
    elif opt in [ '-s', '--statusonly' ]:
        verbose = 0


IN = '''<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/">
<SOAP-ENV:Body>
    <hello/>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>
'''

IN = '''<SOAP-ENV:Envelope
 xmlns="http://www.example.com/schemas/TEST"
 xmlns:SOAP-ENV="http://schemas.xmlsoap.org/soap/envelope/"
 xmlns:xsd="http://www.w3.org/2001/XMLSchema"
 xmlns:SOAP-ENC="http://schemas.xmlsoap.org/soap/encoding/">
<SOAP-ENV:Body>
    <echo>
        <SOAP-ENC:int>1</SOAP-ENC:int>
        <SOAP-ENC:int>2</SOAP-ENC:int>
    </echo>
</SOAP-ENV:Body>
</SOAP-ENV:Envelope>
'''

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((hostname, portnum))
f = s.makefile('r+')

f.write('POST /cgi-bin/x HTTP/1.0\r\n')
f.write('Content-type: text/xml; charset="utf-8"\r\n')
f.write('Content-Length: %d\r\n\r\n' % len(IN))
f.write(IN)
f.flush()

status = f.readline()
print status,
while 1:
    l = f.readline()
    if l == '': break
    if verbose: print l,

f.close()
