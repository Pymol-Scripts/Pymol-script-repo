## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   This is the header for the PDB2PQR Source tree.  All desired 
##          components can be found in the PDB2PQR_SUBDIRS variable.
##
## Author:  Todd Dolinsky
## ###########################################################################
AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
SUBDIRS = @PDB2PQR_SUBDIRS@
distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	-rm -rf autom4te.cache
