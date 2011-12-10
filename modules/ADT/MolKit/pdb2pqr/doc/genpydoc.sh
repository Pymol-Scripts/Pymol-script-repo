#!/bin/bash
#
# Script to generate contents of pydoc subdirectory
# Nathan Baker
# $Id: $

export PYTHONPATH=${PYTHONPATH}:`cd ..; pwd`
echo $PYTHONPATH

# Generate file lists
srclist=`ls ../src/*py | grep -v "__init__"`
srclist=`echo ${srclist} | sed -e "s/\.\.\/src\//src./g" | sed -e "s/\.py//g"`
srclist="pdb2pqr ${srclist}"
exlist="copy string sys re getopt math time xml.sax os"
echo ${srclist}

cd pydoc
cat << EOF > index.html
<html>
	<head>
		<title>PDB2PQR Python Documentation</title>
		<link rel="stylesheet" href="http://agave.wustl.edu/css/baker.css" type="text/css">
	</head>
	<body>
		<h1>PDB2PQR Python Documentation</h1>
		<h2>PDB2PQR modules</h2>
		<ul>
EOF

for src in ${srclist}; do 
	echo ${src}
	pydoc -w ${src}
	echo "		<li><a href="${src}.html">${src}</a></li>" >> index.html
done
cat << EOF >> index.html
		</ul>

		<h2>Other modules</h2>
		<ul>
EOF
for src in ${exlist}; do 
	echo ${src}
	pydoc -w ${src}
	echo "		<li><a href="${src}.html">${src}</a></li>" >> index.html
done
cat << EOF >> index.html
	</ul>
	</body>
</html>
EOF
