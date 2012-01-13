## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Todd Dolinsky
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .c .C .o .so .f
AM_CFLAGS        = $(PY_CFLAGS) $(PY_FLAGS)
CCLD             = $(CC)
AM_FFLAGS        = $(PY_FLAGS)
CLEANFILES       = _propkalib.so propka

noinst_PROGRAMS       = _propkalib.so propka
_propkalib_so_SOURCES   = propkalib.c propka.f propkac.c
_propkalib_so_LINK      = $(F77)
_propkalib_so_LDFLAGS   = -o $@ $(PY_SHARED)
_propkalib_so_LDADD     = $(PY_LDFLAGS)
propka_SOURCES          = propka.f
distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	-rm -rf autom4te.cache
