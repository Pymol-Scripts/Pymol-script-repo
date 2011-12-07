/* PQR cleansing code from Julie Mitchell.  Reformats PQR files to meet PDB standards.

	Julie says:
	
	Here is the code that takes the output from pdb2pqr and reformats
	the charge and radius entries.  It is VERY crude.  All it does is
	shave off the last digit, but for my purposes and probably those
	of the Jmol folks, nine thousandths of an Angstrom or unit charge
	isn't a big deal. 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define ss 120

/* global vars */

int main(int argc, char **argv){

	char helpStr[] = "Usage:  fixpqr  infile.pqr  outfile.pqr  \n";
	FILE *infile, *outfile;
	char inName[ss], outName[ss];
	char inStr[ss], fixStr[ss];
		
	/* read user input */
	if (argc < 3){
		fprintf(stderr,helpStr);
		exit(1);
	} else {	
		strcpy(inName,argv[1]);
		strcpy(outName,argv[2]);
	}
	
	/* open input/output files */
	
	infile = fopen(inName,"r");
	
	if (infile == NULL){
		fprintf(stderr,"ERROR:  could not find file %s\n",inName);
		exit(1);
	}
	
	outfile = fopen(outName,"w");
	
	/* read lines of input file */
	
	while (fgets(inStr,ss,infile) != NULL){
		
		strcpy(fixStr,inStr);
		
		/* check to see if line is an atom record */
		
		if  ( (strncmp (inStr,"ATOM",4) == 0) || (strncmp(inStr,"HETATM",6) == 0) )  {
		
//0123456789012345678901234567890123456789012345678901234567890123456789
//ATOM   2783  HG  CYS   178     -27.703   2.860  33.146  0.2068 0.6000
//ATOM   2783  HG  CYS   178     -27.703   2.860  33.146  0.206 0.600

			strcpy(&fixStr[61],&inStr[62]);
			strcpy(&fixStr[67],&inStr[69]);
			fprintf(outfile,"%s",fixStr);
								    			 
		} else {
			fprintf(outfile,"%s",fixStr);
		}
	}
	
	/* close io files */
	fclose(infile);
	fclose(outfile);
	
return 0;}

