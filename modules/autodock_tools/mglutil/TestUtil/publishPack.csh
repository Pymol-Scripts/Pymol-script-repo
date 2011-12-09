#!/bin/tcsh

# 1) go to the mgl version of the packages
# cd /mgl/python/share/lib/python2.5/site-packages/
cd /usr/local/home/sophiec/MGL

# 2) Tag the current mgl tag with the given new tag
cvs tag -r mgl $1 $2

# 3) Update the mgl version of the given package
cvs co $2

# 4) Move the mgl tag to the latest version of the given packages
cvs tag -F mgl $2

