dnl Some PDB2PQR info

dnl ##########################################################################
dnl File: configure.in
dnl Purpose: configure.in script for PDB2PQR
dnl Note:  For use with autoconf-2.56 and above
dnl Author:  Todd Dolinsky and Nathan Baker
dnl ##########################################################################

dnl ##########################################################################
dnl 0. autoconf initialization
dnl ##########################################################################
AC_INIT([PropKa], [1.0.0], [Jan-Jensen@uiowa.edu])
AC_PREREQ(2.56)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT(${HOME})
AC_CANONICAL_SYSTEM

dnl ##########################################################################
dnl 1. automake initialization (uncomment only if rebuilding Makefile.in's)
dnl ##########################################################################
dnl # THIS IS THE VERSION INFORMATION FOR THE PDB2PQR LIBRARY.  ACCORDING TO 
dnl # THE LIBTOOL MANUAL:
dnl #   This flag accepts an argument of the form `current[:revision[:age]]'.
dnl #   So, passing `-version-info 3:12:1' sets current to 3, revision to 12,
dnl #   and age to 1.  If either revision or age are omitted, they default to 
dnl #   0. Also note that age must be less than or equal to the current 
dnl #   interface number.  Here are a set of rules to help you update your
dnl #   library version information:
dnl #      1. Start with version information of `0:0:0' for each libtool 
dnl #         library.
dnl #      2. Update the version information only immediately before a public
dnl #         release of your software. More frequent updates are unnecessary,
dnl #         and only guarantee that the current interface number gets larger
dnl #         faster.  
dnl #      3. If the library source code has changed at all since the last
dnl #         update, then increment revision (`c:r:a' becomes `c:r+1:a').
dnl #      4. If any interfaces have been added, removed, or changed since the
dnl #         last update, increment current, and set revision to 0.
dnl #      5. If any interfaces have been added since the last public release,
dnl #         then increment age.
dnl #      6. If any interfaces have been removed since the last public
dnl #         release, then set age to 0.
dnl #
dnl #   Never try to set the interface numbers so that they correspond to the
dnl #   release number of your package. This is an abuse that only fosters
dnl #   misunderstanding of the purpose of library versions.
AM_INIT_AUTOMAKE(PropKa, 1.0.0)
PROPKA_VERSION="1:0:0"
AC_SUBST(PROPKA_VERSION)
ACLOCAL=":";
AUTOCONF=":";
AUTOMAKE=":";
AUTOHEADER=":";
AC_SUBST(ACLOCAL)
AC_SUBST(AUTOCONF)
AC_SUBST(AUTOMAKE)
AC_SUBST(AUTOHEADER)
dnl ##########################################################################
dnl 2. checks for programs
dnl ##########################################################################
AC_PATH_PROG(py_path, python)
if test -z "${py_path}"; then
   AC_MSG_ERROR([Unable to find python!])
fi
AC_LANG_C
AC_PROG_CC
AC_PROG_F77
AC_LANG_FORTRAN77
dnl AC_F77_LIBRARY_LDFLAGS 
dnl # -----------------------
dnl # HANDLE OSF WEIRDNESS:
dnl #   This includes the inability to link FORTRAN and C code compiled with
dnl #   debugger flags and a g77 that breaks on some FORTRAN code with the
dnl #   -O2 flag set.
dnl # -----------------------
short_os=`echo ${host_os} | sed 's/osf.*/osf/g'`;
if test -z "`echo ${short_os} | sed 's/osf//g'`"; then
  AC_MSG_RESULT([Disabling debugger compiler flags to avert OSF weirdness])
  if test -z ${FFLAGS-''}; then
      FFLAGS="";
  else
     FFLAGS=`echo ${FFLAGS} | sed 's/-g//g'`;
     if test -z "`echo ${F77-f77} | sed 's/g77//g'`"; then
         FFLAGS=`echo ${FFLAGS} | sed 's/-O2/-O/g'`;
     fi
  fi
  if test -z ${CFLAGS-''}; then
      CFLAGS="";
  else
     CFLAGS=`echo ${CFLAGS} | sed 's/-g//g'`;
  fi
  if test -z ${CPPFLAGS-''}; then
      CPPFLAGS="";
  else
     CPPFLAGS=`echo ${CPPFLAGS} | sed 's/-g//g'`;
  fi
  if test -z ${CXXFLAGS-''}; then
      CXXFLAGS="";
  else
     CXXFLAGS=`echo ${CXXFLAGS} | sed 's/-g//g'`;
  fi
fi

FFLAGS="${FLIBS} ${FFLAGS}"

AC_LANG_C
AC_PROG_CPP

AC_MSG_CHECKING([for python include directory])
changequote(<<, >>)dnl
PY_VERSION=`$py_path -c 'import sys; print sys.version[:3]'`   
changequote([, ])dnl
PY_PREFIX=`$py_path -c 'import sys; print sys.prefix'`
if test -f $PY_PREFIX/include/python$PY_VERSION/Python.h; then    
   PY_CFLAGS="-I$PY_PREFIX/include/python$PY_VERSION"
   AC_MSG_RESULT([$PY_CFLAGS])
else
   AC_MSG_ERROR([Unable to find Python.h!])
fi

dnl # Platform specific linking flags

PY_LDFLAGS=""
PY_SHARED=""
PY_FLAGS=""
case $host_os in
   linux*) 
      PY_SHARED="-shared";;
   irix*) 
      PY_SHARED="-shared"
      PY_LDFLAGS="-lftn -lm";;
   sun|solaris**) 
      PY_SHARED="-G"
      PY_LDFLAGS="-lF77 -lm -lM77 -lsunmath";;
   darwin*)
      PY_LDFLAGS="-Wl,-framework -Wl,Python -bundle";;
esac

case $host_cpu in
   x86_64*)
      PY_FLAGS="-fPIC";;
esac

AC_SUBST(PY_SHARED,[$PY_SHARED])
AC_SUBST(PY_LDFLAGS,[$PY_LDFLAGS])
AC_SUBST(PY_CFLAGS,[$PY_CFLAGS])
AC_SUBST(PY_FLAGS, [$PY_FLAGS])

dnl ##########################################################################
dnl 11. generate the makefiles
dnl ##########################################################################

AC_OUTPUT([
    Makefile
])
