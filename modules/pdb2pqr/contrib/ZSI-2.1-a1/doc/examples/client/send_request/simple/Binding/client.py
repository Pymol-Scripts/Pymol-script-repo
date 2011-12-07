#!/usr/bin/env python

from ZSI import Binding


MESSAGE = "Hello from Python!"

def main():
    binding = Binding(url='http://localhost:8080/server.py')
    print ' Sending: %s' % MESSAGE
    response = binding.echo(MESSAGE)
    print 'Response: %s' % MESSAGE


if __name__ == '__main__':
    main()
