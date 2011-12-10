""" Author:  Alan Ezust
 Version: 2.2
 Date: October 30, 2004
 (relpath.py v1 originally by Cimarron Taylor
 from Oreilly/Activestate Python cookbook 2003)

 helper functions for relative paths.
 This package includes rel2abs() and abs2rel(), 
 based on the perl functions from cpan File::Spec

 Version 2.1 fixes/simplifies pathsplit - uses the string split instead of
 a very inefficient recursive routine.
 Also fixed rel2abs to do a normpath on already absolute paths.
 
 """

import os
import os.path
import re

# matches http:// and ftp:// and mailto://
protocolPattern = re.compile(r'^\w+://')
parent = ".." + os.path.sep # use urlparse.urljoin for relative urls


def isabs(string):
    """ 
    
    @return true if string is an absolute path or protocoladdress
    for addresses beginning in http:// or ftp:// or ldap:// - 
    they are considered "absolute" paths.
    """
    if protocolPattern.match(string): return 1
    return os.path.isabs(string)


def rel2abs(path, base = os.curdir):
    """ converts a relative path to an absolute path.

    @param path the path to convert - if already absolute, is returned
    normalized
    @param base - optional. Defaults to the current location
    The base is intelligently concatenated to the given relative path.
    @return the relative path of path from base
    """
    if isabs(path): return os.path.normpath(path)
    retval = os.path.join(base,path)
    return os.path.abspath(retval)


def commonpath(l1, l2, common=[]):
    if len(l1) < 1: return (common, l1, l2)
    if len(l2) < 1: return (common, l1, l2)
    if l1[0] != l2[0]: return (common, l1, l2)
    return commonpath(l1[1:], l2[1:], common+[l1[0]])


def relpath(base, path):
    """ returns the relative path from base to path """
    baselist = base.split(os.path.sep)
    pathlist = path.split(os.path.sep)
    (common,l1,l2) = commonpath(baselist, pathlist)
    p = []
    if len(l1) > 0:
        p = [ parent * len(l1) ]
    p = p + l2
    if len(p) is 0:
        return "."
    return os.path.join( *p )
    
    
def abs2rel(path, base = os.curdir):
    """ @return a relative path from base to path.
    
    base can be absolute, or relative to curdir, or defaults
    to curdir.
    """
    if protocolPattern.match(path): return path
    base = rel2abs(base)
    return relpath(base, path)


if __name__ == "__main__" : 
    filename = "/home/alan/public_html/oopdocbook/icons/home.png"
    path1 = "/home/alan/public_html/oopdocbook"
    path2 = "/home/alan/public_html/oopdocbook/"

    rel1 = abs2rel(filename, path1)
    rel2 = abs2rel(filename, path2)
    assert (rel1 == rel2)
    
    assert (rel2 == "icons/home.png")
    
    path3 = "/home/alan/public_html/photos/misc/jewel.png"
    path4 = "/home/alan/public_html/oopdocbook/docs/"    
    relpath = abs2rel(path3, path4)
    
    expected = os.path.join("..", "..", "photos", "misc", "jewel.png")
    assert (relpath == expected)
    
    print "done"
    