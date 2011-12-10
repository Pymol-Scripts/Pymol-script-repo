import sys
import string

VERSION = string.split(sys.version)[0]

if VERSION >= '2.1':
    import testplus
else:
    import testplus_before_2_1
    testplus = testplus_before_2_1

