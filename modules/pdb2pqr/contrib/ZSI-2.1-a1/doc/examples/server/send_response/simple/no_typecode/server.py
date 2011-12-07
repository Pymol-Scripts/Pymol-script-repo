#!/usr/bin/env python

from ZSI import dispatch


def echo(message):
    return message


if __name__ == '__main__':
    dispatch.AsServer(port=8080)
