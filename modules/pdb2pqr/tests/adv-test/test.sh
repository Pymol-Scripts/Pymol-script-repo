PY_PATH=/usr/bin/python
SRCDIR=/home/tlinnet/Software/pymol/Pymol-script-repo/modules/pdb2pqr 
CWD=$SRCDIR/tests/adv-test
echo ----------------------------------------------------------------------------------
echo Running pdb2pqr with pH = 7.00, ligand = $CWD/LIG_1HPX.mol2, force field = AMBER, and path = 1HPX
$PY_PATH $SRCDIR/pdb2pqr.py --with-ph=7.00 --ligand=$SRCDIR/tests/adv-test/LIG_1HPX.mol2 --ff=AMBER 1HPX $SRCDIR/tests/adv-test/test-output-user.pqr

$CWD/advtest
