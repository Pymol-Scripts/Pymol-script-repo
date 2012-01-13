#!/usr/bin/env python

from Registration_services import *

import sys

def main():
    user = ns1.User_Def()
    user._UserId = 'john_doe'
    user._Name = 'John Doe'
    user._Age = 25

    locator = RegistrationServiceLocator()
    registration = locator.getRegistration(tracefile=sys.stdout)

    request = RegisterUserRequestWrapper()
    request._User = user
    response = registration.RegisterUser(request)
    print response._Message

if __name__ == '__main__':
    main()
