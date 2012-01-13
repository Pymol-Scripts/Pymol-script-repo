## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   This is the header for the PDB2PQR Source tree.  All desired 
##          components can be found in the PDB2PQR_SUBDIRS variable.
##
## Author:  Todd Dolinsky
## Modified by J. Ren and Yong Huang
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
SUBDIRS = @PDB2PQR_SUBDIRS@
PREFIX = @prefix@
SRCDIR = @srcdir@
PROPKA = @with_propka@
PDB2PKA = @with_pdb2pka@
URL = @URL@
install_dir_var = @install_dir_var@
website_install_var = @website_install_var@
action_install_var = @action_install_var@

CWD=`pwd`
dirpath=$(CWD)/prefixfile 

EXTRA_DIST = *.py *.css dat/* doc/* examples/* html/master-index.html html/README \
	html/server.html.in src/aa.py src/definitions.py src/forcefield.py src/hydrogens.py \
	src/__init__.py src/inputgen.py src/na.py src/pdb.py src/protein.py src/psize.py \
	src/quatfit.py src/routines.py src/server.py.in src/structures.py src/utilities.py \
	tools/* pdb2pka/*.py pdb2pka/*.h pdb2pka/*.i pdb2pka/*.cpp pdb2pka/*.DAT \
	pdb2pka/substruct/*.py pdb2pka/substruct/clique/* pdb2pka/ligandclean \
	extensions/README extensions/template extensions/*.py propka/*.py propka/*.i \
	tests/* pdb2pka/Makefile.in pdb2pka/Makefile.am

distclean-local: clean-local
clean-local: clean-ac
	-rm -rf src/server.py src/aconf.py index.html html/server.html pdb2pqr.cgi tmp/ pdb2pka/*.o pdb2pka/*.so pdb2pka/substruct/*.so pdb2pka/substruct/*.o junk *~ tests/test/test-output-user* tests/test/test.sh tests/adv-test/test-output-user* tests/adv-test/test.sh tests/test-long/test.sh tests/test-webserver/test.sh tests/local-test/test-output-user* tests/local-test/test.sh
	rm -rf querystatus.cgi apbs_cgi.cgi pdb2pqr.py visualize.cgi
	rm -rf *.pyc
	rm -rf tests/adv-test/test.cpp tests/test/test.cpp tests/local-test/test.cpp

clean-ac:
	-rm -rf autom4te.cache

test:
	tests/test/test.sh

adv-test:
	tests/adv-test/test.sh

test-long:
	tests/test-long/test.sh

local-test:
	tests/local-test/test.sh

test-webserver:
	tests/test-webserver/test.sh

install:
	make
	@if test "$(PREFIX)" == "$(CWD)/" ; then \
	  echo Your current source directory is the install directory.  Nothing to be done. ; \
	else \
	  c=`pwd`; \
	  echo $(PREFIX) > prefixfile ; \
	  awk -F"/" '{print NF-1}' ${dirpath} > fieldfile ; \
	  fields=`cat fieldfile`; \
	  awk -F"/" '{print $$NF}' ${dirpath} > lastfieldfile ; \
	  lastfield=`cat lastfieldfile`; \
	  if test -z $$lastfield; then \
	    fields=`expr $$fields - 1` ; \
	  fi ; \
	  dirpath=`pwd`; \
	  cd / ; \
	  awk -F"/" '{print $$2}' $$dirpath/prefixfile > $$dirpath/dir1file ; \
	  dir1=`cat $$dirpath/dir1file` ; \
	  if test -d $$dir1; then \
	    cd $$dir1 ; \
	  else \
	    if test -w .; then \
	      mkdir $$dir1 ; \
	      cd $$dir1 ; \
	    else \
	      echo You do not have permission to write in $$dir1.  Ask your system administrator to make a symbolic link.  Make sure the symbolic link has the same directory name as your PDB2PQR directory. ; \
	      exit 1  ;\
	    fi ; \
	  fi ;  \
	  if test $$fields -eq 1 ; then \
	    dir=$$dir1 ; \
	  fi ; \
	  for (( i = 3 ; i <= $$fields+1; i++ )) \
	    do \
	    dir=`cat $$dirpath/prefixfile | cut -f$$i -d"/"` ; \
	    if test -d $$dir; then \
	      cd $$dir ;\
	    else \
	      if test -w .; then \
	        mkdir $$dir ; \
	        cd $$dir ; \
	      else \
	        echo You do not have permission to write in $$dir.  Ask your system administrator to make a symoblic link.  Make sure the symbolic link has the same directory name as your PDB2PQR directory ;\
	        exit 2 ; \
	      fi ; \
	    fi ; \
	  done ; \
	  cd .. ; \
	  chmod 755 $(CWD)/$$dir ; \
	  if ! test -d $(CWD)/$$dir/tmp ; then \
	    echo "Creating tmp/" ; \
	    mkdir $(CWD)/$$dir/tmp ; \
	  fi ; \
	  chmod 777 $(CWD)/$$dir/tmp ; \
	  echo "Adding server.html" ; \
	  cp $$c/html/server.html $(CWD)/$$dir/index.html ; \
	  echo "Adding python files" ; \
	  cp $$c/*.py* $(CWD)/$$dir/. ; \
	  echo "Adding __init__.py" ; \
	  cp $$c/__init__.py $(CWD)/$$dir/__init__.py ; \
	  echo "Adding pdb2pqr.py" ; \
	  cp $$c/pdb2pqr.py $(CWD)/$$dir/. ; \
	  echo "Adding pdb2pqr.cgi" ; \
	  cp $$c/pdb2pqr.cgi $(CWD)/$$dir/. ; \
	  echo "Adding pdb2pqr.css" ; \
	  cp $$c/apbs_cgi.cgi $(CWD)/$$dir/. ; \
	  echo "Adding apbs_cgi.cgi" ; \
	  cp $$c/visualize.cgi $(CWD)/$$dir/. ; \
	  echo "Adding visualize.cgi" ; \
	  cp $$c/querystatus.cgi $(CWD)/$$dir/. ; \
	  echo "Adding querystatus.cgi" ; \
	  cp $$c/pdb2pqr.css $(CWD)/$$dir/. ; \
	  if ! test -d $(CWD)/$$dir/dat  ; then \
	    mkdir $(CWD)/$$dir/dat ; \
	  fi ; \
	  echo "Adding dat/" ; \
	  cp $$c/dat/*.xml $(CWD)/$$dir/dat/. ; \
	  cp $$c/dat/*.DAT $(CWD)/$$dir/dat/. ; \
	  cp $$c/dat/*.names $(CWD)/$$dir/dat/. ; \
	  if ! test -d $(CWD)/$$dir/doc  ; then \
	    mkdir $(CWD)/$$dir/doc ; \
	  fi ; \
	  echo "Adding doc/" ; \
	  cp $$c/doc/*.html $(CWD)/$$dir/doc/. ; \
	  cp $$c/doc/*.sh $(CWD)/$$dir/doc/. ; \
	  if ! test -d $(CWD)/$$dir/doc/images  ; then \
	    mkdir $(CWD)/$$dir/doc/images ; \
	  fi ; \
	  cp $$c/doc/images/*.png $(CWD)/$$dir/doc/images/. ; \
	  if ! test -d $(CWD)/$$dir/doc/pydoc ; then \
	    mkdir $(CWD)/$$dir/doc/pydoc ; \
	  fi ; \
	  cp $$c/doc/pydoc/*.html $(CWD)/$$dir/doc/pydoc/. ; \
	  if ! test -d $(CWD)/$$dir/examples ; then \
	    mkdir $(CWD)/$$dir/examples ; \
	  fi ; \
	  echo "Adding examples/" ; \
	  cp $$c/examples/*.html $(CWD)/$$dir/examples/. ; \
	  if ! test -d $(CWD)/$$dir/examples/1a1p ; then \
	    mkdir $(CWD)/$$dir/examples/1a1p ; \
	  fi ; \
	  cp $$c/examples/1a1p/*.pdb $(CWD)/$$dir/examples/1a1p/. ; \
	  cp $$c/examples/1a1p/README $(CWD)/$$dir/examples/1a1p/. ; \
	  if ! test -d $(CWD)/$$dir/examples/ligands ; then \
	    mkdir $(CWD)/$$dir/examples/ligands ; \
	  fi ; \
	  cp $$c/examples/ligands/*.mol* $(CWD)/$$dir/examples/ligands/. ; \
	  if ! test -d $(CWD)/$$dir/extensions ; then \
	    mkdir $(CWD)/$$dir/extensions ; \
	  fi ; \
	  echo "Adding extensions" ; \
	  cp $$c/extensions/*.py $(CWD)/$$dir/extensions/. ; \
	  if ! test -d $(CWD)/$$dir/html  ; then \
	    mkdir $(CWD)/$$dir/html ; \
	  fi ; \
	  cp $$c/html/master-index.html $(CWD)/$$dir/html/. ; \
	  cp $$c/html/server.html $(CWD)/$$dir/html/server.html ; \
	  if ! test -d $(CWD)/$$dir/jmol ; then \
	    mkdir $(CWD)/$$dir/jmol ; \
	  fi ; \
	  echo "Adding jmol/" ; \
	  cp -rf $$c/jmol/*.* $(CWD)/$$dir/jmol/. ; \
	  if ! test -d $(CWD)/$$dir/src ; then \
	    mkdir $(CWD)/$$dir/src ; \
	  fi ; \
	  echo "Adding python files from src/" ; \
	  cp $$c/src/*.py* $(CWD)/$$dir/src/. ; \
	  cp $$c/src/server.py $(CWD)/$$dir/src/server.py ; \
	  if test "$(PROPKA)" == "enable" ; then \
	    if ! test -d $(CWD)/$$dir/propka ; then \
	      mkdir $(CWD)/$$dir/propka ; \
	    fi ; \
	    echo "Adding propka"; \
	    cp $$c/propka/propka $(CWD)/$$dir/propka/. ; \
	    cp $$c/propka/propkalib.i $(CWD)/$$dir/propka/. ; \
	    cp $$c/propka/*.py $(CWD)/$$dir/propka/. ; \
	    cp $$c/propka/*.so $(CWD)/$$dir/propka/. ; \
	  fi ; \
	  if test "$(PDB2PKA)" == "enable" ; then \
	    if ! test -d $(CWD)/$$dir/pdb2pka ; then \
	      mkdir $(CWD)/$$dir/pdb2pka ; \
	    fi ; \
	    echo "Adding pdb2pka" ; \
	    cp $$c/pdb2pka/*.so $(CWD)/$$dir/pdb2pka/. ; \
	    cp $$c/pdb2pka/*.py $(CWD)/$$dir/pdb2pka/. ; \
	    cp $$c/pdb2pka/*.DAT $(CWD)/$$dir/pdb2pka/. ; \
	    cp $$c/pdb2pka/*.h $(CWD)/$$dir/pdb2pka/. ; \
	    cp $$c/pdb2pka/*.i $(CWD)/$$dir/pdb2pka/. ; \
	    if ! test -d $(CWD)/$$dir/pdb2pka/ligandclean ; then \
	      mkdir $(CWD)/$$dir/pdb2pka/ligandclean ; \
	    fi ; \
	    cp -rf $$c/pdb2pka/ligandclean/*.py $(CWD)/$$dir/pdb2pka/ligandclean/. ; \
	    if ! test -d $(CWD)/$$dir/pdb2pka/ligandclean/testcomplexes ; then \
	      mkdir $(CWD)/$$dir/pdb2pka/ligandclean/testcomplexes ; \
	    fi ; \
	    cp -rf $$c/pdb2pka/ligandclean/testcomplexes/*.pdb $(CWD)/$$dir/pdb2pka/ligandclean/testcomplexes/. ; \
	    cp -rf $$c/pdb2pka/ligandclean/testcomplexes/*.mol* $(CWD)/$$dir/pdb2pka/ligandclean/testcomplexes/. ; \
	    if ! test -d $(CWD)/$$dir/pdb2pka/substruct ; then \
	      mkdir $(CWD)/$$dir/pdb2pka/substruct ; \
	    fi ; \
	    cp $$c/pdb2pka/substruct/*.py $(CWD)/$$dir/pdb2pka/substruct/. ; \
	    cp $$c/pdb2pka/substruct/*.so $(CWD)/$$dir/pdb2pka/substruct/. ; \
	    if ! test -d $(CWD)/$$dir/pdb2pka/substruct/clique ; then \
	      mkdir $(CWD)/$$dir/pdb2pka/substruct/clique ; \
	    fi ; \
	    cp $$c/pdb2pka/substruct/clique/*.h $(CWD)/$$dir/pdb2pka/substruct/clique/. ; \
	  fi ; \
	  echo PDB2PQR has been installed in $(CWD)/$$dir ; \
	  echo ; \
	  cd $(CWD)/$$dir ; \
	  cd $$c; \
	  rm lastfieldfile prefixfile dir1file fieldfile ; \
	  echo $(URL) > url_file ; \
	  awk -F"/" '{print $$(NF-1)}' url_file > url_file_2 ; \
	  last_url=`cat url_file_2` ; \
	  echo $(PREFIX) > p_file ; \
	  awk -F"/" '{print $$(NF-1)}' p_file > p_file_2 ; \
	  last_p=`cat p_file_2` ; \
	  rm -f url_file* p_file* ; \
	  if test `whoami` = "root" ; then \
            webviewable="true" ; \
	    echo ${prefix} > prefix_file_tmp ; \
            first=`awk -F"/" '{ print $$2 }' prefix_file_tmp` ; \
	    second=`awk -F"/" '{ print $$3 }' prefix_file_tmp` ; \
	    third=`awk -F"/" '{ print $$4 }' prefix_file_tmp` ; \
	    rm -f prefix_file_tmp ; \
	    if test $$first != var || test $$second != www || test $$third != html ; then  \
	      webviewable=false ; \
	    fi ; \
	    if test "$$last_url" != "$$last_p" || test "$$webviewable" = "false" ; then \
	      echo Need to make symbolic link between $(PREFIX) and /var/www/html/"$$last_url" ; \
	      if test -d "/var/www/html/$$last_url" ; then \
	        echo "***" Warning! Directory /var/www/html/"$$last_url" already exists. ; \
		echo Do you want to overwrite this? \(Y/N\); \
		read choice ; \
	        while (test $$choice != "Y" && test $$choice != "y" && test $$choice != "N" && test $$choice != "n")  \
	          do  \
	          echo You have entered $$choice. Please enter Y or N: ; \
	          read choice ; \
	        done ; \
	        if test $$choice = "N" || test $$choice = "n" ; then  \
	          echo Please run configure again and specify a different --prefix ; \
	          exit 0 ; \
	        fi ; \
	  	rm -rf /var/www/html/"$$last_url" ; \
	      fi ; \
	      ln -s $(PREFIX) /var/www/html/"$$last_url" ; \
	      echo "Symbolic link created." ; \
	    fi ; \
	  else  \
	    echo Ask your system administrator to make a symbolic link in /var/www/html by executing the following command; \
	    echo ln -s $(PREFIX) /var/www/html/"$$last_url" ; \
	  fi ; \
	fi


