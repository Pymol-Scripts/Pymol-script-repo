#Create named selections using regular expressions for the protein sequence
 
import pymol
import re
 
aa = { 'ASP' : 'D' , 'GLU' : 'E' , 'GLN' : 'Q' , 'ASN' : 'N' , 'SER' : 'S' ,
       'THR' : 'T' , 'CYS' : 'C' , 'HIS' : 'H' , 'ARG' : 'R' , 'LYS' : 'K' ,
       'MET' : 'M' , 'ALA' : 'A' , 'ILE' : 'I' , 'LEU' : 'L' , 'VAL' : 'V' ,
       'GLY' : 'G' , 'PRO' : 'P' , 'TRP' : 'W' , 'PHE' : 'F' , 'TYR' : 'Y' ,
       'SCY' : 'U' , 'ASX' : 'B' , 'GLX' : 'Z' , 'XXX' : 'X'}
 
#made this before the sequence view option, probably another way to do it now
 
def seqoneint(model):
   pymol.stored.seq = []
   cmd.iterate("%s and name ca"%model,"stored.seq.append(resn)")
   seq = ""
   for x in pymol.stored.seq:
      if aa.has_key(x):
         res = aa[x]
         seq = seq+res
      else:
         seq = seq + '-'
   return seq
 
 
 
def grepsel(model="(all)",stretch="",prefix="",combined="0",single="1"):
   '''
DESCRIPTION
 
    Create selections matching motifs, using python regular expression syntax.
    Motif is automatically converted to uppercase. Motif selections are labelled
    as "prefix_motif_###", where ### is the index for the first residue of the
    match. Prefix defaults to selection name. combined = 1 creates one selection
    for all occurences. single = 1 creates one selection for each occurance
    (the default).
 
USAGE
 
    grepsel selection, motif, [prefix, [combined, [single ]]]
 
EXAMPLES
 
    Create selections for all motifs matching "ESS" (selection_ESS_###,...):
    grepsel selection, ess
 
    Create selections for the PxGY motif with prefix m (m_P.CY_###,...):
    grepsel selection, p.gy, m
    '''
 
   if selection == "(all)":
      selection = "all"
   if prefix == "":
      prefix=selection
 
   stretch = stretch.upper() 
   seq = seqoneint(selection)
   pymol.stored.resi = []
   pymol.stored.chain = []
   cmd.iterate("%s and name ca"%selection,"stored.resi.append(resi);stored.chain.append(chain)")
   motif = re.compile(stretch)
   occurrences = motif.finditer(seq)
   stretchmod = stretch.replace("+","\+")
   stretchmod = stretchmod.replace("?","\?")
 
   print stretchmod
   if combined == "1":
      cmd.select("%s_%s"%(prefix,stretch), "none")
 
 
   for find in occurrences:      
 
      mb = pymol.stored.resi[find.start()]
      me = pymol.stored.resi[find.end()-1]
 
      ch = pymol.stored.chain[find.start()]
      cmd.select("%s_%s_%s%s"%(prefix,stretch,me,ch), "chain %s and (i; %s-%s)"%(ch,int(mb),int(me)))
      if combined == "1":
         cmd.select("%s_%s"%(prefix,stretch),"\"%s_%s\" | (%s and chain %s and (i; %s-%s))"%(prefix,stretchmod,selection,ch,int(mb),int(me)))
 
   cmd.select("none")
   cmd.delete("sel*")
 
cmd.extend("grepsel",grepsel)

