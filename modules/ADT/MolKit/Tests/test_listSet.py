#############################################################################
#
# Author: Sophie I Coon, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2003
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/Tests/test_listSet.py,v 1.2 2003/08/29 17:50:16 sophiec Exp $
#
# $Id: test_listSet.py,v 1.2 2003/08/29 17:50:16 sophiec Exp $
#



import sys
from mglutil.regression import testplus

## class Object:
##     def __init__(self, name='ball', color='blue', shape='spherical'):
##         self.name = name
##         self.color = color
##         self.shape = shape
    

def test_ListSet_Create_dataNoneElementTypeNone():
    from MolKit.listSet import ListSet
    testList = ListSet()
    assert testList.data == []

def test_ListSet_2():
    class Object:
        def __init__(self, name='ball', color='blue', shape='spherical'):
            self.name = name
            self.color = color
            self.shape = shape
    b = Object()
    c = Object(name='cube', color = 'red', shape='square')
    t = Object(name='triangle', color='yellow', shape='triangle')
    d = Object(name='disk', color='silver', shape='round')
    from MolKit.listSet import ListSet
    testList = ListSet([b,c,t], elementType=Object)
    assert len(testList) == 3
    assert testList.color == ['blue', 'red', 'yellow']
    testList.color = ['green', 'yellow', 'cyan', ]
    assert b.color == 'green' and c.color=='yellow' and t.color == 'cyan'

    testList.append(d)
    assert len(testList) == 4
    assert testList.color == ['green', 'yellow', 'cyan', 'silver']
    
def test_ListSet_loopoveremptyset():
    from MolKit.listSet import ListSet
    testList = ListSet([])
    for l in testList:
        print l.__str__()
    del testList

harness = testplus.TestHarness( __name__,
                                funs = testplus.testcollect( globals()),
                                )

if __name__ == '__main__':
    testplus.chdir()
    print harness
    sys.exit( len( harness))
