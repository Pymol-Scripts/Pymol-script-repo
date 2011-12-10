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
AC_INIT([PDB2PQR], [1.2.1], [todd@ccb.wustl.edu])
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
AM_INIT_AUTOMAKE(PDB2PQR, 1.2.1)
PDB2PQR_VERSION="1:2:1"
AC_SUBST(PDB2PQR_VERSION)
ACLOCAL=":";
AUTOCONF=":";
AUTOMAKE=":";
AUTOHEADER=":";
AC_SUBST(ACLOCAL)
AC_SUBST(AUTOCONF)
AC_SUBST(AUTOMAKE)
AC_SUBST(AUTOHEADER)
AC_PROG_CXX

dnl ##########################################################################
dnl 2. checks for programs
dnl ##########################################################################

AC_PATH_PROG(py_path, python)
if test -z "${py_path}"; then
   AC_MSG_ERROR([Unable to find python!])
fi

AC_MSG_CHECKING([for python include directory])
changequote(<<, >>)dnl
PY_VERSION=`$py_path -c 'import sys; print sys.version[:3]'`   
changequote([, ])dnl
PY_PREFIX=`$py_path -c 'import sys; print sys.prefix'`
if test -f $PY_PREFIX/include/python$PY_VERSION/Python.h; then    
   PY_HEADER="-I$PY_PREFIX/include/python$PY_VERSION"
   AC_MSG_RESULT([$PY_HEADER])
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
AC_SUBST(PY_HEADER,[$PY_HEADER])
AC_SUBST(PY_FLAGS, [$PY_FLAGS])

dnl ##########################################################################
dnl 3. Check for available subdirectories
dnl ##########################################################################

PDB2PQR_SUBDIRS=""

# A.  propka

AC_MSG_CHECKING([for --disable-propka])
AC_ARG_ENABLE(propka, 
     [AC_HELP_STRING([--disable-propka],[Disable propka compilation])],
     [ with_propka="disable" ],
     [ with_propka="enable" ]
)
if test "${with_propka}" = "enable"; then
     PDB2PQR_SUBDIRS="propka"
     AC_MSG_RESULT([no])
     AC_MSG_RESULT([])
     AC_MSG_RESULT([------------- Begin Processing PropKa Configure Options -----------])
     cd propka
     ./configure 
     cd ..
     AC_MSG_RESULT([------------- End Processing PropKa Configure Options -----------])
     AC_MSG_RESULT([])
else
     AC_MSG_RESULT([yes])
     PDB2PQR_SUBDIRS=""
fi

# B. pdb2pka

AC_MSG_CHECKING([for --enable-pdb2pka])
AC_ARG_ENABLE(pdb2pka,
     [AC_HELP_STRING([--enable-pdb2pka], [Enable pdb2pka compilation])],
     [ with_pdb2pka="enable" ],
     [ with_pdb2pka="disable" ]
)

if test "$with_pdb2pka" = "enable" ; then
     PDB2PQR_SUBDIRS="$PDB2PQR_SUBDIRS pdb2pka"
     AC_MSG_RESULT([yes])

     AC_MSG_RESULT([])
     AC_MSG_RESULT([------------- Begin Processing pdb2pka Configure Options -----------])

     # Test for Numeric install
     PY_NUMERIC=""
     AC_MSG_CHECKING([for Numeric include files])
     if test -f $PY_PREFIX/include/python$PY_VERSION/Numeric/arrayobject.h; then
        PY_NUMERIC="-I$PY_PREFIX/include/python$PY_VERSION/Numeric"   
        AC_MSG_RESULT([$PY_NUMERIC])
        PY_CFLAGS="$PY_HEADER $PY_NUMERIC"
     fi

     # If failed, test for numpy
     if test "$PY_NUMERIC" = "" ; then
         AC_MSG_CHECKING([for numpy includes files])  
         if test -f $PY_PREFIX/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy/arrayobject.h; then
             PY_NUMERIC="-I$PY_PREFIX/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy"
             AC_MSG_RESULT([$PY_NUMERIC])
             PY_CFLAGS="$PY_HEADER $PY_NUMERIC"
         else
             AC_MSG_ERROR([Unable to find Numeric or numpy installed!])
         fi
     fi

     # Test for APBS Python wrappers
     #AC_MSG_CHECKING([if --with-apbsdir was used])
     #AC_ARG_WITH(apbsdir,
     #   [AC_HELP_STRING([--with-apbsdir=PATH], [full path to top-level APBS directory])],
     #   [apbsdir="$with_apbsdir/tools/python"; AC_MSG_RESULT([yes])],
     #   [apbsdir=""; AC_MSG_RESULT([no])]
     #)

     #if test "$apbsdir" = ""; then
     #   AC_MSG_CHECKING([for APBSDIR environment variable])
     #   if test "$APBSDIR" != ""; then
     #       apbsdir="$APBSDIR/tools/python"
     #       AC_MSG_RESULT([yes])
     #   else
     #       AC_MSG_RESULT([no])
     #       apbsdir="./pdb2pka"
     #   fi
     #fi
  
     #AC_MSG_CHECKING([for compiled APBS Python wrappers])
     #apbserror=0
     #if test -f "$apbsdir/apbslib.py"; then
     #    if test -f "$apbsdir/_apbslib.so"; then
     #        AC_MSG_RESULT([$apbsdir])
     #    else
     #        AC_MSG_RESULT([no])
     #        apbserror=1
     #    fi
     #else
     #    AC_MSG_RESULT([no])
     #    apbserror=1
     #fi

     #if test "$apbserror" -eq 1; then
     #    echo ""
     #    echo "Unable to find apbslib.py and _apbslib.so in $apbsdir!"
     #    echo "Compiled APBS Python wrappers are needed to use the pdb2pka functions"
     #    echo "Please recompile APBS with the '--enable-python' flag."
     #    echo ""
     #    AC_MSG_ERROR(["APBS Python wrappers not found!"])
     #fi


     #AC_SUBST(APBS_PYDIR,[$apbsdir])
     AC_MSG_RESULT([------------- End Processing pdb2ka Configure Options -----------])
     AC_MSG_RESULT([])


else
     PY_CFLAGS="$PY_HEADER"
     AC_MSG_RESULT([no])
fi

AC_SUBST(PY_CFLAGS,[$PY_CFLAGS])
AC_SUBST(PDB2PQR_SUBDIRS,[$PDB2PQR_SUBDIRS])

dnl ##########################################################################
dnl 11. generate the makefiles
dnl ##########################################################################

AC_MSG_RESULT([------------- End Processing PDB2PQR Configure Options -----------])
AC_MSG_RESULT([])
if test -n "${PDB2PQR_SUBDIRS}"; then
   AC_OUTPUT([
        Makefile
        pdb2pka/Makefile
        pdb2pka/substruct/Makefile
   ])
   AC_MSG_RESULT([Please type "make" to complete the installation.])
fi
