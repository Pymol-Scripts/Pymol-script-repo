#!/usr/bin/env python

from ZSI import dispatch

from Registration_services import GetUserResponseWrapper
from Registration_services_types import ns1


def GetUser(user_id):
    user = ns1.User_Def()
    user._UserId = user_id
    user._Name = "John Doe"
    user._Age = 32

    response = GetUserResponseWrapper()
    response._User = user
    return response


if __name__ == '__main__':
    dispatch.AsServer(port=8080)
