#!/usr/bin/env python

import Example_services

MESSAGE = "Hello from Python!"


def main():
    locator = Example_services.ExampleServiceLocator()
    port = locator.getExample()
    request = Example_services.EchoRequestWrapper()
    request._Message = MESSAGE
    
    try:
        print ' Sending: %s' % MESSAGE
        response = port.echo(request)
        print 'Response: %s' % response._Message
    except Exception, e:
        print e


if __name__ == '__main__':
    main()
