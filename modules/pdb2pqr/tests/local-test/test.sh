PY_PATH=/usr/bin/python
SRCDIR=/home/tlinnet/Software/pymol/Pymol-script-repo/modules/pdb2pqr 
echo ----------------------------------------------------------------------------------
echo This is a local command line test.  Running pdb2pqr with forcefield=AMBER and path=1AFS.pdb
$PY_PATH $SRCDIR/pdb2pqr.py --ff=AMBER $SRCDIR/tests/local-test/1AFS.pdb $SRCDIR/tests/local-test/test-output-user.pqr

CWD=$SRCDIR/tests/local-test
$CWD/localtest

