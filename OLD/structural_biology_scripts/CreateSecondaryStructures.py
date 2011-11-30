##############################################
#   Original Author:  Dan Kulp
#   Date  :  9/8/2005
#    MOdified by Jurgen F. Doreleijers
#    For Hamid Eghbalnia               
#
#############################################
# Call in window like : 
# @C:\Documents and Settings\jurgen.WHELK.000\workspace\Wattos\python\Wattos\Utils\CreateSecondaryStructures.py
# Next line is a pymol directive
python
import string
import urllib
 
 
# Well I guess one can build a protein with it but the vdw contacts would be horrible.
# Peptide needs to be at least 2 residues.
def createPeptide(seqInfo):
    cmd.delete("all")
    # Creates residue TWO
    editor.attach_amino_acid('pk1',seqInfo[1][0]) 
    # Creates residue ONE
    createSS('resi 2', sequence=seqInfo[0][0],terminal='N')
    print "found sequence info for number of residues: ", len(seqInfo)
    for i in range(2,len(seqInfo) ):
        # resn is the residue number of the new residue
        resn = i + 1
        print "Adding residue: ", resn,   seqInfo[i][0]
        # Note that the previous residue is numbered i. 
        resi = 'resi '+`i`
        createSS(resi, sequence=seqInfo[i][0])
    for i in range( len(seqInfo) ):
        resi = 'resi '+`i+1`
#        print "Setting backbone angles for residue: ", (i+1),   seqInfo[i][0],seqInfo[i][1],seqInfo[i][2]
        set_phipsi(resi,seqInfo[i][1],seqInfo[i][2])
 
# Create generic secondary structure, based off a selection
def createSS(sel, sequence='ALA',repeat=1,terminal='C'):
 
    # Set selection
    selection = "%s and name %s" % (sel,terminal)
 
    # Pick atom for editing - interestingly only need to do this for the first addition
    cmd.edit(selection,None,None,None,pkresi=0,pkbond=0)
 
    # Array of residues
    seq = string.split(sequence,",")
 
    # Get residue numbering .. potential bug here if number is inconsistent.. (Only works at c-terminal)
    resi = int(cmd.get_model(sel).atom[0].resi) + 1
    # Loop and build new residues
    for i in range(1,repeat+1):
        for s in seq:
#            print "residue[%i]: %s %s" % (i,s,terminal)
            editor.attach_amino_acid('pk1',s)
 
    # Remove extra OXT carboxylate atom (OXT1, OXT2 ?) .. fix as needed
    if terminal == 'C':
        cmd.remove("%s and name OXT" % sel)
 
 
def set_phipsi(sel,phi,psi):
    # Get atoms from selection
    atoms = cmd.get_model("byres ("+sel+")")
 
    # Loop through atoms in selection        
    for at in atoms.atom:
        if at.name == "N":
            # Check for a null chain id (some PDBs contain this) 
            unit_select = ""
            if not at.chain == "":
               unit_select = "chain "+str(at.chain)+" and "
 
            try:
                # Define residue selections     
                residue_def_prev = unit_select+'resi '+str(int(at.resi)-1)
                residue_def      = unit_select+'resi '+str(at.resi)        
#                print "residue_def_prev: [%s]" % residue_def_prev
#                print "residue_def     : [%s]" % residue_def
                if at.resn == "PRO":
                    print "Skipping setting phi for PRO"
                else:
                    old_phi = cmd.get_dihedral(residue_def_prev+' and name C',residue_def+' and name N', residue_def+' and name CA',residue_def+' and name C')        
                    cmd.set_dihedral(          residue_def_prev+' and name C',residue_def+' and name N', residue_def+' and name CA',residue_def+' and name C'      ,phi)
                    print "Changed residue %4s %4s phi: from %6.1f to %6.1f" % (at.resn, at.resi, old_phi, float(phi))        
            except:
 
                print "Note skipping set of phi because of error; this is normal for a N-terminal residue"
            try:
                residue_def      = unit_select+'resi '+str(at.resi)
                residue_def_next = unit_select+'resi '+str(int(at.resi)+1)
#                print "residue_def     : [%s]" % residue_def
#                print "residue_def_next: [%s]" % residue_def_next
                old_psi = cmd.get_dihedral(residue_def     +' and name N',residue_def+' and name CA',residue_def+' and name C', residue_def_next+' and name N')
                cmd.set_dihedral(          residue_def     +' and name N',residue_def+' and name CA',residue_def+' and name C', residue_def_next+' and name N',psi)
                print "Changed residue %4s %4s psi: from %6.1f to %6.1f" % (at.resn, at.resi, old_psi, float(psi))        
            except:
                print "Note skipping set of psi; this is normal for a C terminal residue"
 
def getTableFromCsvFile(urlLocation):
  result = []
  r1 = urllib.urlopen(urlLocation)
  data = r1.read()
  r1.close()  
  dataLines = data.split("\n")   
  for dataLine in dataLines:
    if dataLine:
        result.append( dataLine.split(',') )     
  return result
 
# next line is a pymol directive.
python end
 
os.chdir("C:\Documents and Settings\jurgen.WHELK.000\workspace\Wattos\python\Wattos\Utils")
seqInfo = getTableFromCsvFile("seqInfo.csv")
createPeptide(seqInfo)

