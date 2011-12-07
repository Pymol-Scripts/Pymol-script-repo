#!/usr/bin/env python

from ZSI import dispatch

from Example_services import EchoResponseWrapper


def echo(message):
    response = EchoResponseWrapper()
    response._Message = message
    return response


if __name__ == '__main__':
    dispatch.AsServer(port=8080)
