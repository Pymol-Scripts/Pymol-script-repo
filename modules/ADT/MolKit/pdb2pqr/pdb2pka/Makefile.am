## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   This is the driver for the pdb2pka directory.
##
## Author:  Todd Dolinsky
## ###########################################################################
AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .H .h .C .cc .c .f .o .cpp .so
AM_CXXFLAGS        = $(PY_CFLAGS) $(PY_FLAGS)
CCLD             = $(CXX)
CLEANFILES       = _pMC_mult.so *.pyc apbslib.py _apbslib.so
SUBDIRS          = substruct

noinst_PROGRAMS         = _pMC_mult.so
_pMC_mult_so_SOURCES     = pMC_mult.cpp pMC_mult_wrap.cpp  
_pMC_mult_so_LINK        = $(CXX)
_pMC_mult_so_LDFLAGS     = -o $@ $(PY_SHARED)
_pMC_mult_so_LDADD       = $(PY_LDFLAGS)

#all:
#	ln -s -f $(APBS_PYDIR)/apbslib.py .
#	ln -s -f $(APBS_PYDIR)/_apbslib.so .

distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	-rm -rf autom4te.cache
