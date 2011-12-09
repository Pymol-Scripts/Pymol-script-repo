## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   This is the driver for the substruct directory of pdb2pka.
##
## Author:  Todd Dolinsky
## ###########################################################################
AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .H .h .C .cc .c .f .o .cpp .so
AM_CXXFLAGS        = $(PY_CFLAGS) $(PY_FLAGS)
CCLD             = $(CXX)
CLEANFILES       = Algorithms.so 

noinst_PROGRAMS          = Algorithms.so
Algorithms_so_SOURCES     = Algorithms.cpp  
Algorithms_so_LINK        = $(CXX)
Algorithms_so_LDFLAGS     = -o $@ $(PY_SHARED)
Algorithms_so_LDADD       = $(PY_LDFLAGS)

distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	-rm -rf autom4te.cache
