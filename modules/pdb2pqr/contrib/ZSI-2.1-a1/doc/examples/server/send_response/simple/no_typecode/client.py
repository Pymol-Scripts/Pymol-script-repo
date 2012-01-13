#!/usr/bin/env python

from ZSI import ServiceProxy
import sys

MESSAGE = "Hello from Python!"

def main():
    server = ServiceProxy('../binding.wsdl', use_wsdl=False)

    print ' Sending: %s' % MESSAGE
    response = server.echo(MESSAGE)
    print 'Response: %s' % response


if __name__ == '__main__':
    main()
