#! /usr/bin/env python
# $Header$
'''Simple CGI dispatching.
'''

from ZSI import *
from ZSI import _copyright
import base64, os

_b64_decode = base64.decodestring

# Typecode to parse a ZSI BasicAuth header.
_auth_tc = TC.Struct(None,
                        [ TC.String('Name'), TC.String('Password') ],
                        extras=1)

class AUTH:
    '''Constants for authentication mechanisms.
    '''
    none = 0
    httpbasic = 1
    zsibasic = 2
    httpdigest = 4

class ClientBinding:
    '''Information about the client that is connected to us.
    '''

    def __init__(self, ps):
        self.ps, self.auth = \
            ps, None
        self.environ = os.environ.copy()
        self.environ['CONTENT_LENGTH'] = str(0)

    def GetAuth(self):
        '''Return a tuple containing client authentication data.
        '''
        if self.auth: return self.auth
        for elt in self.ps.GetMyHeaderElements():
            if elt.localName == 'BasicAuth' \
            and elt.namespaceURI == ZSI_SCHEMA_URI:
                d = _auth_tc.parse(elt, self.ps)
                self.auth = (AUTH.zsibasic, d['Name'], d['Password'])
                return self.auth
        ba = self.environ.get('HTTP_AUTHENTICATION')
        if ba:
            ba = ba.split(' ')
            if len(ba) == 2 and ba[0].lower() == 'basic':
                ba = _b64_decode(ba[1])
                self.auth = (AUTH.httpbasic,) + tuple(ba.split(':'))
                return self.auth
        self.auth = (AUTH.none,)
        return self.auth

    def GetNS(self):
        '''Return namespace for the top main request element.
        '''
        return self.ps.body_root.namespaceURI or ''

    def GetRequest(self):
        '''Return the ParsedSoap request.
        '''
        return self.ps

if __name__ == '__main__': print _copyright
