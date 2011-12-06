#!/usr/bin/env python

from ZSI import dispatch

from ComplexTypes import User

def GetUser(user_id):
    user = User(user_id, 'John Doe', 28)
    return user

if __name__ == '__main__':
    nsdict = { 'types' : 'http://pycon.org/types' }
    dispatch.AsServer(port=8080, nsdict=nsdict)
