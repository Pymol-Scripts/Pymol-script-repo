#!/usr/bin/env python

from ZSI import dispatch
import Registration_services_types
from Registration_services import RegisterUserResponseWrapper

def RegisterUser(user):
    response = RegisterUserResponseWrapper()
    response._Message = "OK"
    return response

if __name__ == '__main__':
    dispatch.AsServer(port=8080, typesmodule=(Registration_services_types,))
