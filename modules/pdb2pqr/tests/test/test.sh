PY_PATH=/usr/bin/python
SRCDIR=/home/tlinnet/Software/pymol/Pymol-script-repo/modules/pdb2pqr 
echo ----------------------------------------------------------------------------------
echo This is a simple command line test.  Running pdb2pqr with forcefield=AMBER and path=1AFS
$PY_PATH $SRCDIR/pdb2pqr.py --ff=AMBER 1AFS $SRCDIR/tests/test/test-output-user.pqr

CWD=$SRCDIR/tests/test
$CWD/simpletest

echo You may also run "make adv-test" to test propka and pdb2pka.
