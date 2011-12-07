#!/usr/bin/env python

from ZSI import dispatch

import ComplexTypes as MyComplexTypes


def RegisterUser(user):
    return "OK"


if __name__ == '__main__':
    dispatch.AsServer(port=8080, typesmodule=(MyComplexTypes,))
