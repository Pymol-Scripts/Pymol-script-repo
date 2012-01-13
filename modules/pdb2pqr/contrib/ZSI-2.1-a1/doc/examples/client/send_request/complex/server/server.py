#!/usr/bin/env python

from ZSI import dispatch

import ComplexTypes as MyComplexTypes

from Registration_services import RegisterUserResponseWrapper


def RegisterUser(user):
    response = RegisterUserResponseWrapper()
    response._Message = "OK"
    return response


if __name__ == '__main__':
    dispatch.AsServer(port=8080, typesmodule=(MyComplexTypes,))
