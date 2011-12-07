#! /usr/bin/env python
# $Header$
'''Utilities for HTTP Digest Authentication
'''
import re
from md5 import md5
import random
import time
import httplib

random.seed(int(time.time()*10))

def H(val):
  return md5(val).hexdigest()

def KD(secret,data):
  return H('%s:%s' % (secret,data))

def A1(username,realm,passwd,nonce=None,cnonce=None):
  if nonce and cnonce:
    return '%s:%s:%s:%s:%s' % (username,realm,passwd,nonce,cnonce)
  else:
    return '%s:%s:%s' % (username,realm,passwd)

def A2(method,uri):
  return '%s:%s' % (method,uri)

def dict_fetch(d,k,defval=None):
  if d.has_key(k):
    return d[k]
  return defval

def generate_response(chaldict,uri,username,passwd,method='GET',cnonce=None):
  """
  Generate an authorization response dictionary. chaldict should contain the digest
  challenge in dict form. Use fetch_challenge to create a chaldict from a HTTPResponse
  object like this: fetch_challenge(res.getheaders()).

  returns dict (the authdict)

  Note. Use build_authorization_arg() to turn an authdict into the final Authorization
  header value.
  """
  authdict = {} 
  qop = dict_fetch(chaldict,'qop')
  domain = dict_fetch(chaldict,'domain')
  nonce = dict_fetch(chaldict,'nonce')
  stale = dict_fetch(chaldict,'stale')
  algorithm = dict_fetch(chaldict,'algorithm','MD5')
  realm = dict_fetch(chaldict,'realm','MD5')
  opaque = dict_fetch(chaldict,'opaque')
  nc = "00000001"
  if not cnonce:
    cnonce = H(str(random.randint(0,10000000)))[:16]

  if algorithm.lower()=='md5-sess':
    a1 = A1(username,realm,passwd,nonce,cnonce)
  else:
    a1 = A1(username,realm,passwd)

  a2 = A2(method,uri)

  secret = H(a1)
  data = '%s:%s:%s:%s:%s' % (nonce,nc,cnonce,qop,H(a2))
  authdict['username'] = '"%s"' % username
  authdict['realm'] = '"%s"' % realm
  authdict['nonce'] = '"%s"' % nonce
  authdict['uri'] = '"%s"' % uri
  authdict['response'] = '"%s"' % KD(secret,data)
  authdict['qop'] = '"%s"' % qop
  authdict['nc'] = nc
  authdict['cnonce'] = '"%s"' % cnonce
  
  return authdict


def fetch_challenge(http_header):
  """ apparently keywords Basic and Digest are not being checked 
  anywhere and decisions are being made based on authorization 
  configuration of client, so I guess you better know what you are
  doing.  Here I am requiring one or the other be specified.

      challenge Basic auth_param
      challenge Digest auth_param
  """
  m = fetch_challenge.wwwauth_header_re.match(http_header)
  if m is None:
      raise RuntimeError, 'expecting "WWW-Authenticate header [Basic,Digest]"'

  d = dict(challenge=m.groups()[0])
  m = fetch_challenge.auth_param_re.search(http_header)
  while m is not None:
      k,v = http_header[m.start():m.end()].split('=')
      d[k.lower()] = v[1:-1]
      m = fetch_challenge.auth_param_re.search(http_header, m.end())
       
  return d

fetch_challenge.wwwauth_header_re = re.compile(r'\s*([bB]asic|[dD]igest)\s+(?:[\w]+="[^"]+",?\s*)?')
fetch_challenge.auth_param_re = re.compile(r'[\w]+="[^"]+"')


def build_authorization_arg(authdict):
  """
  Create an "Authorization" header value from an authdict (created by generate_response()).
  """
  vallist = []
  for k in authdict.keys():
    vallist += ['%s=%s' % (k,authdict[k])]
  return 'Digest '+', '.join(vallist)

if __name__ == '__main__': print _copyright
