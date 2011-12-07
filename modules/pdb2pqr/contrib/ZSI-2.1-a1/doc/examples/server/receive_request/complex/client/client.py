#!/usr/bin/env python

from ZSI import ServiceProxy

from ComplexTypes import User

import sys

def main():
    user = User('john_doe', 'John Doe', 25)
    nsdict = { 'types' : 'http://pycon.org/typs' }
    registration = ServiceProxy('../binding.wsdl',
                                nsdict=nsdict,
                                tracefile=sys.stdout)
    response = registration.RegisterUser(user)
    print response

if __name__ == '__main__':
    main()
