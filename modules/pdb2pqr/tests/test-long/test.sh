PY_PATH=/usr/bin/python
SRCDIR=/home/tlinnet/Software/pymol/Pymol-script-repo/modules/pdb2pqr 
echo "----------------------------------------------------------------------------------"
echo "This is a command line test.  Running PDB2PQR on a long list of PDBs (246 PDBs)."
echo "However, if you prefer to run a shorter list of PDBs, e.g., 50 PDBs, you can try"
echo "the following commands instead:"
echo "  export TESTNUM=50; make test-long" 
echo "----------------------------------------------------------------------------------"

if test "$TESTNUM" = ""; then
  $PY_PATH $SRCDIR/tests/test-long/test.py
else
  $PY_PATH $SRCDIR/tests/test-long/test.py --testnum=`echo $TESTNUM`
fi
