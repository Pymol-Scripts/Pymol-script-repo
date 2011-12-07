#!/usr/bin/env python

from ZSI import ServiceProxy

MESSAGE = "Hello from Python!"

def main():
    server = ServiceProxy('../binding.wsdl', use_wsdl=True)

    print ' Sending: %s' % MESSAGE
    response = server.echo(Message=MESSAGE)
    print 'Response: %s' % response['Message']


if __name__ == '__main__':
    main()
