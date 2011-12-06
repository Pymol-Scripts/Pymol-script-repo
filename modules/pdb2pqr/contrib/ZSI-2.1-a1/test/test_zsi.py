#!/usr/bin/env python
import unittest
import test_t1
import test_t2
import test_t3
import test_t5
import test_t6
import test_t7
import test_t8
import test_t9
import test_union
import test_list
import test_TCtimes
import test_URI
import test_rfc2617
import test_QName

def makeTestSuite():
    return unittest.TestSuite(
        map(lambda t: globals()[t].makeTestSuite(), 
            filter(lambda g: g.startswith('test_') and True, globals()))
    )

def main():
    unittest.main(defaultTest="makeTestSuite")
    suite = unittest.TestSuite()

if __name__ == "__main__" : main()
