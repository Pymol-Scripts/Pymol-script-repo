#!/usr/bin/env python

import sys

from Registration_services import *

def main():
    locator = RegistrationServiceLocator()
    service = locator.getRegistration(tracefile=sys.stdout)

    request = GetUserRequestWrapper()
    request._UserId = "john_doe"

    response = service.GetUser(request)
    print '   Age: %d' % response._User._Age
    print '  Name: %s' % response._User._Name
    print 'UserId: %s' % response._User._UserId
    

if __name__ == '__main__':
    main()
