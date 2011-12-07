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
AC_INIT([PDB2PQR], [1.7], [baker@biochem.wustl.edu])
AC_PREREQ(2.56)
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_SYSTEM

# Added by J.Ren

chmod 755 .

pwd > cwdfile
cwd=`awk -F"/" '{ print $NF }' cwdfile`
rm cwdfile

if test ${prefix} == "NONE" ; then
  if test `whoami` = "root" ; then
    if test -d "/var/www/html/pdb2pqr" ; then
      echo "*** Warning: As root, your default --prefix is set to /var/www/html/pdb2pqr."
      echo "    However, we have detected that /var/www/html/pdb2pqr already exists."
      echo "    You may run configure again and specify a different --prefix"
      echo "    Do you want to overwrite this directory? (Y/N)"
      read choice 
      while test $choice != "Y" && test $choice != "y" && test $choice != "N" && test $choice != "n"
        do
        echo You have entered $choice. Please enter Y or N:
	read choice
      done 	
      if test $choice = "N" || test $choice = "n" ; then
        echo Please run configure again and specify a different --prefix
	exit 0
      fi  
    fi
    prefix=/var/www/html/pdb2pqr
  else
    if test -d ${HOME}/pdb2pqr ; then
      echo "*** Warning: As a regular user, your default --prefix is set to ${HOME}/pdb2pqr."
      echo "    However, we have detected that ${HOME}/pdb2pqr already exists."
      echo "    You may run configure again and specify --prefix"
      echo "    You may run configure again and specify a different --prefix"
      echo "    Do you want to overwrite this directory? (Y/N)"
      read choice 
      while test $choice != "Y" && test $choice != "y" && test $choice != "N" && test $choice != "n"
        do
        echo You have entered $choice. Please enter Y or N:
	read choice
      done 	
      if test $choice = "N" || test $choice = "n" ; then
        echo Please run configure again and specify a different --prefix
	exit 0
      fi  
    fi
    prefix=${HOME}/pdb2pqr
  fi
  default_prefix_used="true"
fi

echo ${prefix} > prefix_test 
isslash=`awk -F"/" '{print $NF}' prefix_test`
prefix=`sed -e 's|//|/|g' prefix_test`
rm -f prefix_test

if ! test -z $isslash ; then 
  prefix=${prefix}/
fi

AC_SUBST(prefix)
AC_PROG_INSTALL

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
AM_INIT_AUTOMAKE(PDB2PQR, 1.7)
PDB2PQR_VERSION="1:7"
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


AC_ARG_WITH(python,
  [  --with-python=PYTHONPATH           Specify Python location.],
  [PYTHONPATH=$withval],
  [PYTHONPATH=""])

if test "$PYTHONPATH" = "" ; then
  AC_PATH_PROG(py_path, python) 
  if test -z "${py_path}"; then
     AC_MSG_ERROR([Unable to find python!  Python location can be specified with --with-python])
  fi
else
  py_path=$PYTHONPATH	
fi

srcpath=`pwd`

sed -e 's%@srcdir@%'${srcpath}'%' tests/test/test.cpp.in > tests/test/test.cpp
sed -e 's%@srcdir@%'${srcpath}'%' tests/adv-test/test.cpp.in > tests/adv-test/test.cpp
sed -e 's%@srcdir@%'${srcpath}'%' tests/local-test/test.cpp.in > tests/local-test/test.cpp
sed -e 's%@py_path@%'${py_path}'%' tests/test/test.sh.in > tests/test/test.sh.tmp.1
sed -e 's%@py_path@%'${py_path}'%' tests/adv-test/test.sh.in > tests/adv-test/test.sh.tmp.1
sed -e 's%@py_path@%'${py_path}'%' tests/test-long/test.sh.in > tests/test-long/test.sh.tmp.1
sed -e 's%@py_path@%'${py_path}'%' tests/test-webserver/test.sh.in > tests/test-webserver/test.sh.tmp.1
sed -e 's%@py_path@%'${py_path}'%' tests/local-test/test.sh.in > tests/local-test/test.sh.tmp.1
sed -e 's%@srcdir@%'${srcpath}'%' tests/test/test.sh.tmp.1 > tests/test/test.sh
sed -e 's%@srcdir@%'${srcpath}'%' tests/adv-test/test.sh.tmp.1 > tests/adv-test/test.sh
sed -e 's%@srcdir@%'${srcpath}'%' tests/test-long/test.sh.tmp.1 > tests/test-long/test.sh
sed -e 's%@srcdir@%'${srcpath}'%' tests/test-webserver/test.sh.tmp.1 > tests/test-webserver/test.sh
sed -e 's%@srcdir@%'${srcpath}'%' tests/local-test/test.sh.tmp.1 > tests/local-test/test.sh
rm -f tests/test/test.sh.tmp.*
rm -f tests/adv-test/test.sh.tmp.*
rm -f tests/test-long/test.sh.tmp.*
rm -f tests/test-webserver/test.sh.tmp.*
rm -f tests/local-test/test.sh.tmp.*

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
     [AC_HELP_STRING([--disable-propka], [Disable propka compilation])],
     [ with_propka="$enableval" ],
     [ with_propka="enable" ]
)

if test "$with_propka" = "yes" ; then
  with_propka="enable"
elif test "$with_propka" = "no"; then
  with_propka="disable"
fi

if test "${with_propka}" = "enable"; then
     PDB2PQR_SUBDIRS="propka"
     AC_MSG_RESULT([no])
     AC_MSG_RESULT([])
     AC_MSG_RESULT([------------- Begin Processing PropKa Configure Options -----------])
     cd propka
     ./configure py_path=$py_path
     cd ..
     AC_MSG_RESULT([------------- End Processing PropKa Configure Options -----------])
     AC_MSG_RESULT([])
else
     AC_MSG_RESULT([yes])
     PDB2PQR_SUBDIRS=""
fi

# B. pdb2pka

AC_MSG_CHECKING([for--disable-pdb2pka])
AC_ARG_ENABLE(pdb2pka,
     [AC_HELP_STRING([--disable-pdb2pka], [Disable pdb2pka compilation])],
     [ with_pdb2pka="$enableval" ],
     [ with_pdb2pka="enable" ]
)

if test "$with_pdb2pka" = "yes" ; then
  with_pdb2pka="enable"
elif test "$with_pdb2pka" = "no"; then
  with_pdb2pka="disable"
fi

if test "$with_pdb2pka" = "enable" ; then
  PY_CFLAGS="$PY_HEADER"
  PACKAGE_PATH=""
  PDB2PQR_SUBDIRS="$PDB2PQR_SUBDIRS pdb2pka"
  AC_MSG_RESULT([no])

  AC_MSG_RESULT([])
  AC_MSG_RESULT([------------- Begin Processing pdb2pka Configure Options -----------])

dnl #######################
dnl  Modified by Yong Huang
dnl #######################

  # Test for Numeric install
  AC_MSG_CHECKING([whether you defined Numeric or Numpy])
  if (test -n "${NUMERIC}") ; then
    AC_MSG_RESULT([yes])
    NUMERIC_INC=`find ${NUMERIC} -iname "arrayobject.h" | tail -n1`
    if test "$NUMERIC_INC" = ""; then
      AC_MSG_ERROR([Please check your NUMERIC path.  arrayobject.h cannot be found in ${NUMERIC}])
    fi
    NUMERIC_INC=`dirname ${NUMERIC_INC}`
    NUMERIC_LIB=`find ${NUMERIC} -iname "Numeric.pth" | tail -n1`
    if test "${NUMERIC_LIB}" = ""; then
      AC_MSG_ERROR([Please check your NUMERIC path.  Numeric.pth cannot be found in ${NUMERIC}])
    fi
    NUMERIC_NAME=`cat ${NUMERIC_LIB}`
    NUMERIC_LIB=`dirname ${NUMERIC_LIB}`
    NUMERIC_LIB=${NUMERIC_LIB}/${NUMERIC_NAME}
    if test "${PACKAGE_PATH}" = ""; then
      PACKAGE_PATH=${NUMERIC_LIB}
    else
      PACKAGE_PATH="${PACKAGE_PATH}:${NUMERIC_LIB}"
    fi
    PY_NUMERIC="-I${NUMERIC_INC} "
    PY_CFLAGS="$PY_CFLAGS $PY_NUMERIC"
  else
    AC_MSG_RESULT([no])	
    AC_MSG_CHECKING([for Numeric include files])
    if test -f $PY_PREFIX/include/python$PY_VERSION/Numeric/arrayobject.h; then
      PY_NUMERIC="-I$PY_PREFIX/include/python$PY_VERSION/Numeric"   
      AC_MSG_RESULT([$PY_NUMERIC])
      PY_CFLAGS="$PY_CFLAGS $PY_NUMERIC"
    else
      AC_MSG_RESULT([no])
    fi
  fi

  # If failed, test for numpy
  AC_MSG_CHECKING([whether you defined numpy])
  if (test -n "${NUMPY}") ; then
    AC_MSG_RESULT([yes])
    NUMPY_INC=`find ${NUMPY} -iname "arrayobject.h" | tail -n1`
    if test "$NUMPY_INC" = ""; then
        AC_MSG_ERROR([Please check your NUMPY path.  arrayobject.h cannot be found in ${NUMPY}])
    fi
    NUMPY_INC=`dirname ${NUMPY_INC}`
    NUMPY_LIB=`find ${NUMPY} -iname "site-packages" | tail -n1`
    if test "${NUMPY_LIB}" = ""; then
        AC_MSG_ERROR([Please check your NUMPY path.  site-packages cannot be found in ${NUMPY}])
    fi
    NUMPY_LIB=`dirname ${NUMPY_LIB}`
    NUMPY_LIB=${NUMPY_LIB}/site-packages
    if test "${PACKAGE_PATH}" = ""; then
      PACKAGE_PATH=${NUMPY_LIB}
    else
      PACKAGE_PATH="${PACKAGE_PATH}:${NUMPY_LIB}"
    fi
    PY_NUMPY="-I${NUMPY_INC}"
    PY_CFLAGS="$PY_CFLAGS $PY_NUMPY"
  else
    AC_MSG_RESULT([no]) 
    if test "$PY_NUMERIC" = "" ; then
      AC_MSG_CHECKING([for numpy includes files])  
      if test -f $PY_PREFIX/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy/arrayobject.h; then
        PY_NUMERIC="-I$PY_PREFIX/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy"
        AC_MSG_RESULT([$PY_NUMERIC])
        PY_CFLAGS="$PY_CFLAGS $PY_NUMERIC"
      elif test -f $PY_PREFIX/Extras/lib/python/numpy/core/include/numpy/arrayobject.h; then
        PY_NUMERIC="-I$PY_PREFIX/Extras/lib/python/numpy/core/include/numpy"
        AC_MSG_RESULT([$PY_NUMERIC])
        PY_CFLAGS="$PY_CFLAGS $PY_NUMERIC"
      elif test -f ${prefix}numpy/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy/arrayobject.h; then
        PY_NUMERIC="-I${prefix}numpy/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy"
        AC_MSG_RESULT([$PY_NUMERIC])
        PY_CFLAGS="$PY_CFLAGS $PY_NUMERIC"
        SITEPKGS=${prefix}numpy/lib/python$PY_VERSION/site-packages
        if test "${PACKAGE_PATH}" = ""; then
          PACKAGE_PATH=${SITEPKGS}
        else
          PACKAGE_PATH="${PACKAGE_PATH}:${SITEPKGS}"
        fi
      elif test -f /opt/rocks/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy/arrayobject.h; then
        AC_MSG_RESULT([$PY_NUMERIC])
        PY_FLAGS="$PY_HEADER -I/opt/rocks/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy"
      else
        $py_path ./contrib/numpy-1.1.0/setup.py install --prefix=${prefix}numpy
        PY_NUMERIC="-I${prefix}numpy/lib/python$PY_VERSION/site-packages/numpy/core/include/numpy"
        AC_MSG_RESULT([$PY_NUMERIC])
        PY_CFLAGS="$PY_CFLAGS $PY_NUMERIC"
        SITEPKGS=${prefix}numpy/lib/python$PY_VERSION/site-packages
        if test "${PACKAGE_PATH}" = ""; then
          PACKAGE_PATH=${SITEPKGS}
        else
          PACKAGE_PATH="${PACKAGE_PATH}:${SITEPKGS}"
        fi
      fi
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
  AC_MSG_RESULT([yes])
fi

#PY_CFLAGS="$PY_CFLAGS -I$PY_PREFIX/lib/python$PY_VERSION "

AC_SUBST(PY_CFLAGS,[$PY_CFLAGS])
AC_SUBST(PDB2PQR_SUBDIRS,[$PDB2PQR_SUBDIRS])

dnl #################
dnl  Added by J Ren
dnl #################

AC_ARG_WITH(url,
  [  --with-url=URL           Set url for the website.  Default http://fully_qualified_domain/pdb2pqr.],
  [URL=$withval],
  [URL="http://`uname -n`/pdb2pqr"])

echo ${URL} > url_test 
isslash=`awk -F"/" '{print $NF}' url_test`
hashttp=`awk -F"http://" '{print $2}' url_test`
rm -f url_test

if ! test -z $isslash ; then 
  URL=${URL}/
fi

if test -z $hashttp ; then
  URL=http://${URL}
fi

AC_SUBST(URL)

dnl ####################
dnl   Added by W Goodman
dnl ####################

dnl ########################
dnl   Modified by Yong Huang
dnl ########################

AC_ARG_WITH(max-atoms,
  [  --with-max-atoms=MAX_ATOMS    Sets the maximum number of atoms in a protein for non-Opal job submission.  Default is 10000],
  [MAX_ATOMS=$withval],
  [MAX_ATOMS=10000])

HAVE_PDB2PQR_OPAL=0

AC_ARG_WITH(opal,
	[  --with-opal=PDB2PQR_OPAL_URL			Set url for the Opal service.  Default is http://kryptonite.nbcr.net/opal2/services/pdb2pqr_1.7.1],
	[PDB2PQR_OPAL_URL=$withval],
	[PDB2PQR_OPAL_URL=""])

if test "$PDB2PQR_OPAL_URL" = "yes"; then
  PDB2PQR_OPAL_URL="http://kryptonite.nbcr.net/opal2/services/pdb2pqr_1.7.1"
  HAVE_PDB2PQR_OPAL=1
fi

if test "$PDB2PQR_OPAL_URL" != ""; then 
  HAVE_PDB2PQR_OPAL=1
  if test -d $PY_PREFIX/lib/python$PY_VERSION/site-packages/ZSI ; then
    echo ZSI is installed
  else
    #  need to install ZSI
    echo Installing ZSI.  This might take a while...
    cd contrib/ZSI-2.1-a1 
    if test `whoami` = "root"; then
      $py_path setup.py install
    else
      zsi_install="true"
      mkdir -p ${prefix}lib/python$PY_VERSION/site-packages
      export PYTHONPATH=${prefix}lib/python$PY_VERSION/site-packages
      $py_path setup.py install --prefix=${prefix}
      if test -d ${prefix}ZSI ; then
        echo ZSI is installed
      else
        cp -r ./build/lib/ZSI ${prefix}ZSI
      fi
      PYTHONPATH=""
    fi 
    cd ../..
  fi
fi  

HAVE_APBS=""
APBS_OPAL_URL=""

AC_ARG_WITH(apbs,
	[  --with-apbs=APBSPATH       Specify APBS executable location.],
	[APBSPATH=$withval],
	[APBSPATH=""])

if test "$APBSPATH" != ""; then
  HAVE_APBS=$APBSPATH
fi

# Passing apbslib.py and _apbslib.so files from APBS to pdb2pka
echo ${HAVE_APBS} > HAVE_APBS_FILE 
HAVE_APBS_PYTHON_DIR=`sed -e 's|bin/apbs|share/tools/python|g' HAVE_APBS_FILE`
HAVE_APBS_MESH_DIR=`sed -e 's|bin/apbs|share/tools/mesh|g' HAVE_APBS_FILE`
if test -e ${HAVE_APBS_PYTHON_DIR}/apbslib.py; then
  cp  ${HAVE_APBS_PYTHON_DIR}/apbslib.py ./pdb2pka
fi
if test -e ${HAVE_APBS_PYTHON_DIR}/_apbslib.so; then
  cp ${HAVE_APBS_PYTHON_DIR}/_apbslib.so ./pdb2pka
fi
if test -e ${HAVE_APBS_MESH_DIR}/smooth; then
  cp  ${HAVE_APBS_MESH_DIR}/smooth ./pdb2pka
fi
rm -f HAVE_APBS_FILE

AC_ARG_WITH(apbs-opal,
	[  --with-apbs-opal=APBS_OPAL_URL       Specify APBS Opal url. Default is http://kryptonite.nbcr.net/opal2/services/apbs_1.3],
	[APBS_OPAL_URL=$withval],
	[APBS_OPAL_URL=""])

if test "$APBS_OPAL_URL" = "yes"; then
  APBS_OPAL_URL="http://kryptonite.nbcr.net/opal2/services/apbs_1.3"
  DEFAULT_APBS_OPAL_URL="1"
else
  DEFAULT_APBS_OPAL_URL="0"
fi


if test "$PDB2PQR_OPAL_URL" = "" && test "$APBS_OPAL_URL" != ""; then
  if test -d $PY_PREFIX/lib/python$PY_VERSION/site-packages/ZSI ; then
    echo ZSI is installed
  else
    #  need to install ZSI
    echo Installing ZSI.  This might take a while...
    cd contrib/ZSI-2.1-a1 
    if test `whoami` = "root"; then
      $py_path setup.py install
    else
      zsi_install="true"
      mkdir -p ${prefix}lib/python$PY_VERSION/site-packages
      export PYTHONPATH=${prefix}lib/python$PY_VERSION/site-packages
      $py_path setup.py install --prefix=${prefix}
      if test -d ${prefix}ZSI ; then
        echo ZSI is installed
      else
        cp -r ./build/lib/ZSI ${prefix}ZSI
      fi
      PYTHONPATH=""
    fi 
    cd ../..
  fi
fi  

# Add tmp/ if it doesn't already exist

AC_CHECK_FILE(tmp/, hastmp="yes", hastmp="no")

if test "${hastmp}" = "no"; then
  mkdir tmp
fi

cp .htaccess tmp

chmod 777 tmp

chmod o+x tests/test/test.sh
chmod u+x tests/test/test.sh
chmod o+x tests/adv-test/test.sh
chmod u+x tests/adv-test/test.sh
chmod o+x tests/test-long/test.sh
chmod u+x tests/test-long/test.sh
chmod o+x tests/test-webserver/test.sh
chmod u+x tests/test-webserver/test.sh
chmod o+x tests/local-test/test.sh
chmod u+x tests/local-test/test.sh

if test "${prefix}" != `pwd` ; then 
  website_install=${URL}
  pdb2pqr_cgi_install=${URL}pdb2pqr.cgi

  echo ${website_install} > website_install_file
  echo ${prefix} > prefix_file
  last_website_install=`awk -F"/" '{print $(NF-1)}' website_install_file`
  last_prefix=`awk -F"/" '{print $(NF-1)}' prefix_file`
  rm -f website_install_file prefix_file
fi

#if test "${prefix}" != `pwd` ; then 
#  if test "${default_prefix_used}" = "true" ; then
#    website_install=${URL}
#    pdb2pqr_cgi_install=${URL}pdb2pqr.cgi
#  else 
#    echo ${prefix} > install_dir_file
#    lastdir=`awk -F"/" '{print $NF}' install_dir_file`

#    if test -z $lastdir; then
#      last=`awk -F"/" '{print NF-1}' install_dir_file`
#      lastdir=`awk -F"/" '{print $(NF-1)}' install_dir_file`
#    fi

#    rm -f install_dir_file

#    website_install=${URL}$lastdir
#    pdb2pqr_cgi_install=${URL}$lastdir/pdb2pqr.cgi
#  fi
#fi

dnl #####################
dnl   Added by Yong Huang
dnl #####################

# Create pdb2pqr.py, main.py and main.cgi; handle the shebangs
sed -e 's%@WHICHPYTHON@%'${py_path}'%' pdb2pqr.py.in > pdb2pqr.py
chmod 755 pdb2pqr.py
cp pdb2pqr.py pdb2pqr.cgi
chmod -x apbs_cgi.py

sed -e 's%@WHICHPYTHON@%'${py_path}'%' apbs_cgi.py > apbs_cgi.cgi
chmod 755 apbs_cgi.cgi

sed -e 's%@WHICHPYTHON@%'${py_path}'%' visualize.py > visualize.cgi
chmod 755 visualize.cgi

#echo ${srcpath} > srcpathtemp1
#sed 's/\//\\\//g' srcpathtemp1 > srcpathtemp2
#srcpathvar=`cat srcpathtemp2`
#rm srcpathtemp1 srcpathtemp2
#sed -e 's/\@srcpath\@/'${srcpathvar}'/g' src/server.py.in > src/server.py.in.tmp.1

#if test "${prefix}" != `pwd` ; then 
  #echo ${prefix} > install_dir_temp1
  #sed 's/\//\\\//g' install_dir_temp1 > install_dir_temp2
  #install_dir_var=`cat install_dir_temp2`
  #rm install_dir_temp1 install_dir_temp2
  #  sed -e 's/\@srcpath\@/'${install_dir_var}'/g' src/server.py.in > src/server.py.install.tmp.1
#fi

echo ${prefix} > install_dir_temp1
sed 's/\//\\\//g' install_dir_temp1 > install_dir_temp2
install_dir_var=`cat install_dir_temp2`
rm install_dir_temp1 install_dir_temp2
sed -e 's/\@srcpath\@/'${install_dir_var}'/g' src/server.py.in > src/server.py.in.tmp.1

AC_SUBST(install_dir_var)

website=${URL}
echo $website > websitetemp1
sed 's/\//\\\//g' websitetemp1 > websitetemp2
websitevar=`cat websitetemp2`
rm websitetemp1 websitetemp2 
sed -e 's/\@website\@/'${websitevar}'/g' src/server.py.in.tmp.1 > src/server.py.in.tmp.2	
mv src/server.py.in.tmp.2 src/server.py
rm -rf src/server.py.in.tmp.* 

#if test "${prefix}" != `pwd` ; then 
#  echo ${website_install} > website_install_temp1
#  sed 's/\//\\\//g' website_install_temp1 > website_install_temp2
#  website_install_var=`cat website_install_temp2`
#  rm website_install_temp1 website_install_temp2 
  # sed -e 's/\@website\@/'${website_install_var}'/g' src/server.py.install.tmp.1 > src/server.py.install.tmp.2	
  # cp src/server.py.install.tmp.2 src/server.py.install
  # rm -rf src/server.py.install.tmp.* 
#fi

AC_SUBST(website_install_var)

pdb2pqr_cgi=${URL}pdb2pqr.cgi
echo $pdb2pqr_cgi > actiontemp1
sed 's/\//\\\//g' actiontemp1 > actiontemp2
actionvar=`cat actiontemp2`
rm actiontemp1 actiontemp2 
if test -z "${PDB2PQR_OPAL_URL}" ; then
  sed -e 's/\@action\@/'${actionvar}'/g' \
      -e '/@WITHOPAL@/d' \
      -e 's%@MAXATOMS@%'${MAX_ATOMS}'%' \
      -e '/@WITHOUTOPAL@/,/@WITHOUTOPAL@/d' html/server.html.in > html/server.html.tmp
  sed -e 's/\@website\@/'${websitevar}'/g' html/server.html.tmp > html/server.html.tmp.2
  cp html/server.html.tmp.2 html/server.html
  rm html/server.html.tmp*
else
  sed -e 's/\@action\@/'${actionvar}'/g' \
      -e '/@WITHOPAL@/,/@WITHOPAL@/d' \
      -e 's/\@website\@/'${websitevar}'/g' \
      -e '/@WITHOUTOPAL@/d' html/server.html.in > html/server.html
fi
cp html/server.html index.html

if test "${prefix}" != `pwd` ; then 
  pdb2pqr_cgi_install=${URL}pdb2pqr.cgi
  echo $pdb2pqr_cgi_install > action_install_temp1
  sed 's/\//\\\//g' action_install_temp1 > action_install_temp2
  action_install_var=`cat action_install_temp2`
  rm action_install_temp1 action_install_temp2 
  #  sed -e 's/\@action\@/'${action_install_var}'/g' html/server.html.in > html/server.html.install.tmp
  #  sed -e 's/\@website\@/'${website_install_var}'/g' html/server.html.install.tmp > html/server.html.install.tmp.2
  #  cp html/server.html.install.tmp.2 html/server.html.install
  #  rm html/server.html.install.tmp*
fi

AC_SUBST(action_install_var)

sed -e 's%@PDB2PQR_OPAL_URL@%'${PDB2PQR_OPAL_URL}'%' \
    -e 's%@WHICHPYTHON@%'${py_path}'%' \
    -e 's%@HAVE_PDB2PQR_OPAL@%'${HAVE_PDB2PQR_OPAL}'%' querystatus.py > querystatus.cgi
sed -e 's%@MAXATOMS@%'${MAX_ATOMS}'%' \
    -e 's%@PACKAGE_PATH@%'${PACKAGE_PATH}'%' \
    -e 's%@PDB2PQR_OPAL_URL@%'${PDB2PQR_OPAL_URL}'%' \
    -e 's%@HAVE_PDB2PQR_OPAL@%'${HAVE_PDB2PQR_OPAL}'%' \
    -e 's%@HAVE_APBS@%'${HAVE_APBS}'%' \
    -e 's%@APBS_OPAL_URL@%'${APBS_OPAL_URL}'%' \
    -e 's%@DEFAULT_APBS_OPAL_URL@%'${DEFAULT_APBS_OPAL_URL}'%' \
    -e 's%@website@%'${website}'%' \
    -e 's%@srcpath@%'${srcpath}'%' \
    -e 's%@INSTALLDIR@%'${prefix}'%' src/aconf.py.in > src/aconf.py
#chmod 755 pdb2pqr-opal.cgi
chmod 755 querystatus.cgi

AC_SUBST(with_propka)
AC_SUBST(with_pdb2pka)

PDB2PQR_SUBDIRS="${PDB2PQR_SUBDIRS} tests"

# Checking for OpenBabel
echo
echo "checking for OpenBabel..."
r=`which babel`
if test ! -z $r; then
  echo "OpenBabel found in $r, you can use OpenBabel to generate MOL2 files."
else
  echo "OpenBabel not found."
fi
echo

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
	tests/Makefile
	tests/test/Makefile
	tests/adv-test/Makefile
	tests/test-long/Makefile
	tests/test-webserver/Makefile
	tests/local-test/Makefile
   ])
   echo
   echo ---------------------------------------------------------------------------------------
   AC_MSG_RESULT([Please run "make" followed by "make install" to complete the installation. You may also run "make test" or "make adv-test". ])
   echo
fi

dnl ##################################
dnl Print out configuration parameters 
dnl ##################################

USER=`whoami`

echo "     ========================"
echo "     Configuration Parameters"
echo "     ========================"
echo
echo "User ID:       ${USER}"
echo "Version:       ${VERSION}"
echo "Install Dir:   ${prefix}"
echo "propka:        ${with_propka}"
echo "pdb2pka:       ${with_pdb2pka}"
echo "PY_CFLAGS:     ${PY_CFLAGS}"
echo 
#echo "     ---------------"
#echo "     Based on SRCDIR"
#echo "     ---------------"
#echo 
#echo "Note:  You cannot access the server if SRCPATH is not a web viewable directory"
#echo
#echo "SRCPATH (location of your source directory):       ${srcpath}"
#echo "WEBSITE (path to the website directory):           ${website}"
#echo "PDB2PQR_CGI (path to pdb2pqr.cgi):                 ${pdb2pqr_cgi}"
#echo
#echo "     ------------------------------"
#echo "     Based on --preifx / INSTALLDIR"
#echo "     ------------------------------"
#echo
#echo "Note:  You cannot access the server if INSTALLDIR is not a web viewable directory"
#echo

if test $prefix = `pwd` ; then
  echo "WEBSITE (path to the website directory):           ${website}"
  echo "   (WEBSITE can be changed with the --with-url option)"
  echo "PDB2PQR_CGI (path to pdb2pqr.cgi):                 ${pdb2pqr_cgi}"
  echo "PDB2PQR_OPAL_URL (url to the pdb2pqr opal service):  ${PDB2PQR_OPAL_URL}"
else
  echo "WEBSITE (path to the website directory):           ${website_install}"
  echo "   (WEBSITE can be changed with the --with-url option)"
  echo "PDB2PQR_CGI (path to pdb2pqr.cgi):                 ${pdb2pqr_cgi_install}"
  if test "$PDB2PQR_OPAL_URL" = ""; then
    echo "PDB2PQR jobs will be forked on the server"
    echo "MAX_ATOMS (max number of atoms in a protein):     ${MAX_ATOMS}"
  else
    echo "PDB2PQR_OPAL_URL (url to the pdb2pqr opal service): ${PDB2PQR_OPAL_URL}"
  fi
  echo
fi
AC_SUBST(last_website_install,[$last_website_install])
AC_SUBST(last_prefix,[$last_prefix])

echo ${prefix} > prefix_file
first=`awk -F"/" '{ print $2 }' prefix_file`
second=`awk -F"/" '{ print $3 }' prefix_file`
third=`awk -F"/" '{ print $4 }' prefix_file`
rm -f prefix_file

if test `whoami` = "root" ; then
  if test ${first} != "var" || test ${second} != "www" || test ${third} != "html" ; then
    echo
    echo "*** Warning: The install directory is not inside /var/www/html/"
    echo "             A symbolic link will be created when you run \"make install\""
  elif test "${last_website_install}" != "${last_prefix}" ; then
    echo
    echo "*** Warning: The last directory of your prefix = ${last_prefix}"
    echo "             The last directory of your URL    = ${last_website_install}"
    echo "             It's recommended that they are the same.  You may run configure again with a different set of --prefix and --with-url"
    echo "             Or  a symbolic link will be created when you run \"make install\", and you will be able to run the server from the URL you specified, but not the URL based on --prefix"
  fi 
else
  echo As a regular user, you must ask your system administrator to create a symbolic link that enables you to view $URL after you run \"make install\"
fi

if test "$zsi_install" = "true" ; then
  echo Since you are not root, ZSI have been installed in your installation directory
fi
