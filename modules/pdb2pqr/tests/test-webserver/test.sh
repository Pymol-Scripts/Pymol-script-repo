PY_PATH=/usr/bin/python
SRCDIR=/home/tlinnet/Software/pymol/Pymol-script-repo/modules/pdb2pqr 
echo "----------------------------------------------------------------------------------"
echo "\nThis is a test for the PDB2PQR web server.  Testing a long list of PDBs (246 PDBs)"
echo "on the PDB2PQR production server.\n"
echo "----------------------------------------------------------------------------------"

$PY_PATH $SRCDIR/tests/test-webserver/web_server_tester.py
