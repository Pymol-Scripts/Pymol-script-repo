'''
Described at PyMOL wiki:
http://www.pymolwiki.org/index.php/aAKMT_Lys_pred

Authors : Troels Schwarz-Linnet
Program : aKMT_Lys_pred
Date    : June 2017

aKMT_Lys_pred -- Help predicting lysine methylation
'''
# Internal pymol import
from pymol import cmd
from pymol import stored
# From Pymol-script repo: https://pymolwiki.org/index.php/Findseq
import findseq

aa_1_3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS',
            'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL',
            'W': 'TRP', 'Y': 'TYR'}

aa_3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

aa_types = {'A': 'hydrophobic', 'C': 'cysteine', 'D': 'negative', 'E': 'negative', 'F': 'aromatic', 'G': 'glycine', 'H': 'polar',
            'I': 'hydrophobic', 'K': 'positive', 'L': 'hydrophobic', 'M': 'hydrophobic', 'N': 'polar', 'P': 'proline', 'Q': 'polar',
            'R': 'positive', 'S': 'polar', 'T': 'polar', 'V': 'hydrophobic', 'W': 'aromatic', 'Y': 'aromatic'}


def get_resis_from_resn(target_sel="all", resn="lys", atom_name="CA", verb=True):
    # Make uppercase
    resn = resn.upper()
    # Check if one letter residue name
    if len(resn) == 1:
        resn = aa_1_3[resn]
    atom_name = atom_name.upper()

    # Prepare for storing and make expression
    stored.infolist = []
    # resv (int): the residue identifier (residue number), ss (str): secondary structure, name (str): the atom name
    expression = "stored.infolist.append([model, chain, resv, resn, ss, name])"
    # Iterate over selection, storing info
    cmd.iterate(target_sel, expression)

    # Store info
    return_list_resn_resi = []
    return_list_resn_resi_sel = []
    group = "Find_%s_%s"%(resn, target_sel)
    for info in stored.infolist:
        cur_model, cur_chain, cur_resv, cur_resn, cur_ss, cur_name = info
        if cur_resn == resn and cur_name == atom_name:
            # Convert residue name to one letter
            cur_aa_3_1 = aa_3_1[resn]
            # Do selection and group
            #sel_str = "/%s/%s//%s"%(cur_model, cur_chain, cur_resv)
            sel_str = "%s and chain %s and resi %s"%(cur_model, cur_chain, cur_resv)
            resn_resi = "%s%s"%(cur_aa_3_1, cur_resv)
            sel_str_text = "%s_%s"%(group, resn_resi)
            cmd.select(sel_str_text, sel_str)
            # Store
            return_list_resn_resi.append(resn_resi)
            return_list_resn_resi_sel.append(sel_str)
            # If verbose
            if verb:
                print("%s , sel: %s"%(resn_resi, sel_str))
    # Group selections
    cmd.group(group, "%s_*"%group)
    cmd.select("%s_sel"%group, "%s_*"%group)
    cmd.show("lines", group)

    # If verbose
    if verb:
        print("\nThere are %i hits, in target_sel=%s, with resn=%s\n"%(len(return_list_resn_resi), target_sel, resn))
    return return_list_resn_resi, return_list_resn_resi_sel
cmd.extend("get_resis_from_resn", get_resis_from_resn)


def match_peptides(target_sel="all", peptides=[], verb=True):
    if type(peptides) != list:
        print("\nERROR: The peptides should be supplied as a list\n")
        return

    # Store info
    return_list_resn_resi = []
    return_list_resn_resi_sel = []
    group = "Match_%s"%(target_sel)

    for peptide in peptides:
        sequence, modification = peptide
        sequence = sequence.strip()
        modification = modification.strip()
        # Check input
        if modification[0].isdigit() or not modification[1:].isdigit():
            print("\nERROR: The modificaions should be in format of ex: K10\n")
            return
        #sel_str = "findseq"
        sel_str = sequence
        findseq.findseq(needle=sequence, haystack=target_sel, selName=sel_str)
        # Limit the selection to atom name CA
        atom_name = "CA"
        sel_str_atom = "%s_%s"%(sel_str, atom_name)
        # Make a sub selection with atom name, and delete old selection
        cmd.select(sel_str_atom, "%s and name %s"%(sel_str, atom_name))
        cmd.delete(sel_str)
        # Iterate
        stored.infolist = []
        # resv (int): the residue identifier (residue number), ss (str): secondary structure, name (str): the atom name
        expression = "stored.infolist.append([model, chain, resv, resn, ss, name])"
        # Iterate over selection, storing info
        cmd.iterate(sel_str_atom, expression)
        cmd.delete(sel_str_atom)

        # Check for results. Is there any found?
        if len(stored.infolist) == 0:
            print("\n#####################################################################")
            print("ERROR: The following sequence cannot be found: %s     %s"%(sequence, modification))
            print("\n#####################################################################\n")
            continue

        # Find resn and index for search. This is for the peptide
        resn = modification[0].upper()
        index = int(modification[1:]) - 1
        # This is for the mathc selection
        peptide_match_modification = stored.infolist[index]
        peptide_match_modification_model, peptide_match_modification_chain, peptide_match_modification_resv, peptide_match_modification_resn, \
            peptide_match_modification_ss, peptide_match_modification_name =  peptide_match_modification
        # Convert to single aa
        peptide_match_modification_resn = aa_3_1[peptide_match_modification_resn]
        # Convert ss, secondary structure, if ss=S (Sheet), or ss='' (Not Helix or Sheet)
        if peptide_match_modification_ss == '':
            peptide_match_modification_ss = 'L'

        # Check if the residue type match
        if peptide_match_modification_resn != resn:
            print("\n#####################################################################")
            print("ERROR: The match is not equal: %s=%s != %s"%(modification[0], resn, peptide_match_modification_resn))
            print("\n#####################################################################\n")
            continue

        # Do selection and group
        #peptide_match_modification_sel_str = "/%s/%s//%s"%(peptide_match_modification_model, peptide_match_modification_chain, peptide_match_modification_resv)
        peptide_match_modification_sel_str = "%s and chain %s and resi %s"%(peptide_match_modification_model, peptide_match_modification_chain, peptide_match_modification_resv)
        peptide_match_modification_resn_resi = "%s%s"%(peptide_match_modification_resn, peptide_match_modification_resv)
        peptide_match_modification_sel_str_text = "%s_%s_%s_%s"%(group, peptide_match_modification_resn_resi, modification, sequence)
        cmd.select(peptide_match_modification_sel_str_text, peptide_match_modification_sel_str)
        # Store
        if peptide_match_modification_resn_resi not in return_list_resn_resi:
            return_list_resn_resi.append(peptide_match_modification_resn_resi)
            return_list_resn_resi_sel.append(peptide_match_modification_sel_str)

        # Print
        if verb:
            print("The peptide=%s, with modification=%s, corresponds to resi=%s"%(sequence, modification, peptide_match_modification_resn_resi))

    # Group selections
    cmd.group(group, "%s_*"%group)
    cmd.select("%s_sel"%group, "%s_*"%group)
    cmd.show("lines", group)

    # If verbose
    if verb:
        print("\nThere are %i uniq matches, in target_sel=%s\n"%(len(return_list_resn_resi), target_sel))
    return return_list_resn_resi, return_list_resn_resi_sel
cmd.extend("match_peptides", match_peptides)


def get_resi_stats(target_sel="all", residues=[], group_id="X", atom="NZ", atom_dist=8, resi_n_term=8, resi_c_term=8,  verb=True):
    # The distance in angstrom to look for
    var_dist = 12

    if type(residues) != list:
        print("\nERROR: The residues should be supplied as a list\n")
        return

    # Get current setting
    ini_setting = cmd.get("dot_solvent")
    ini_setting2 = cmd.get("dot_density")
    # Increasing dot_density makes calculation slower, but not a big difference
    #cmd.set('dot_density', 3)

    # Make groups
    group = "Stats_%s_%s" % (target_sel, group_id)
    group_atom = "Stats_%s_%s_%s" % (atom, target_sel, group_id)
    group_chain = "Stats_%s_%s_%s" % ("chain", target_sel, group_id)
    group_3dweb = "Stats_%s_%s_%s" % ("3dweb", target_sel, group_id)

    # Make list for storing
    slist = []

    # Make file for writing
    wfileweblogo = open("resi_stats_weblogo_%s_%s.txt" % (target_sel, group_id), 'w')

    for residue in residues:
        residue = residue.strip()
        resn_1 = residue[0].upper()
        resn_3 = aa_1_3[resn_1]
        resi = int(residue[1:])
        # Check input
        if resn_1.isdigit():
            print("\nERROR: The residue should be in format of ex: K10\n")
            return

        # Do selection and group
        sel_str = "%s and resn %s and resi %s" % (target_sel, resn_3, resi)
        resn_resi = "%s%s" % (resn_1, resi)
        sel_str_text = "%s_%s_%s" % (target_sel, group_id, resn_resi)
        cmd.select(sel_str_text, sel_str)

        # Make quick test, to see if the atom is there
        sel_str_atom_test = "%s and name %s" % (sel_str_text, atom)
        test_str = "Test_nr_atoms"
        cmd.select(test_str, sel_str_atom_test)
        nr_test = cmd.count_atoms(test_str)
        if nr_test != 1:
            print("\nERROR: The selection '%s', has only nr of atoms:%s. SKIPPING"%(sel_str_atom_test, nr_test))
            continue

        # MSA = Molecular Surface Area
        cmd.set("dot_solvent", "off")
        MSA = cmd.get_area(sel_str)
        # SASA = Solvent Accessible Surface Area
        cmd.set("dot_solvent", "on")
        SASA = cmd.get_area(sel_str)

        # Get the chain residues
        chain = "."*(resi_n_term + resi_c_term + 1)
        chain_sec = "."*(resi_n_term + resi_c_term + 1)
        resi_sel_min = resi-resi_n_term
        if resi_sel_min < 1:
            resi_sel_min = 1
        resi_sel_max = resi+resi_c_term
        resi_sel = "%i-%i" % (resi_sel_min, resi_sel_max)

        # Make selection
        sel_str_chain = "%s and resi %s and name CA" % (target_sel, resi_sel)
        sel_str_text_chain = "%s_%s_%s_%s" % ("chain", target_sel, group_id, resn_resi)
        cmd.select(sel_str_text_chain, sel_str_chain)
        # Get the chain info
        stored.list_chain = []
        expression_chain="stored.list_chain.append([resi, resn, name, ss])"
        cmd.iterate(sel_str_text_chain, expression_chain)
        for chain_resi_info in stored.list_chain:
            chain_resi, chain_resn, chain_name, chain_ss = chain_resi_info
            # Convert ss, secondary structure, if ss=S (Sheet), or ss='' (Not Helix or Sheet)
            if chain_ss == '':
                chain_ss = 'L'
            chain_resi = int(chain_resi)
            try:
                chain_resn_1 = aa_3_1[chain_resn]
            except KeyError:
                chain_resn_1 = "."
            # Calculate index
            index = resi_n_term - (resi - chain_resi)
            # Replace in string for residue names
            chain = chain[:index] + chain_resn_1 + chain[index + 1:]
            # Replace in string for secondary structyre
            chain_sec = chain_sec[:index] + chain_ss + chain_sec[index + 1:]

        # Get number of neighbour atoms
        # Make selection for NZ atoms
        sel_str_atom = "%s and name %s" % (sel_str_text, atom)
        sel_str_text_atom = "%s_%s_%s_%s" % (atom, target_sel, group_id, resn_resi)
        cmd.select(sel_str_text_atom, sel_str_atom)

        # Make selection around NZ atom for fixed distance, and count
        sel_str_atom_around = "%s around %s and not (%s)" % (sel_str_text_atom, atom_dist, sel_str)
        sel_str_text_atom_around = "%s_around_%s_%s_%s" % (atom, target_sel, group_id, resn_resi)
        cmd.select(sel_str_text_atom_around, sel_str_atom_around)
        # Count around
        stored.list = []
        expression="stored.list.append([resi, resn, name])"
        cmd.iterate(sel_str_text_atom_around, expression)
        nr_atoms_around = len(stored.list)

        # Make selection around NZ atom for variable distance
        #for i in range(2, var_dist+1):
        for i in range(2, var_dist+1, 2):
            dist = i
            dist_pre = dist - 1
            # Select for an angstrom shorter
            sel_str_atom_3dweb_pre = "byres %s around %s" % (sel_str_text_atom, dist_pre)
            sel_str_text_atom_3dweb_pre = "%s_3dweb_pre_%s_%s_%s_%s_%s" % (atom, target_sel, group_id, resn_resi, dist, dist_pre)
            cmd.select(sel_str_text_atom_3dweb_pre, sel_str_atom_3dweb_pre)
            # Select at distance
            sel_str_atom_3dweb_post = "byres %s around %s" % (sel_str_text_atom, dist)
            sel_str_text_atom_3dweb_post = "%s_3dweb_post_%s_%s_%s_%s_%s" % (atom, target_sel, group_id, resn_resi, dist, dist)
            cmd.select(sel_str_text_atom_3dweb_post, sel_str_atom_3dweb_post)
            # Make selection for uniq residues with shell
            sel_str_text_atom_3dweb_sel = "%s_3dweb_sel_%s_%s_%s_%s" % (atom, target_sel, group_id, resn_resi, dist)
            cmd.select(sel_str_text_atom_3dweb_sel, "(%s and not %s) and name CA" % (sel_str_atom_3dweb_post, sel_str_atom_3dweb_pre))
            # delete
            cmd.delete(sel_str_text_atom_3dweb_pre)
            cmd.delete(sel_str_text_atom_3dweb_post)
            # Loop through selecion
            stored.list_3dweb = []
            expression_3dweb="stored.list_3dweb.append([resi, resn, name])"
            cmd.iterate(sel_str_text_atom_3dweb_sel, expression_3dweb)
            for web3d_residues in stored.list_3dweb:
                web3d_resi, web3d_resn, web3d_name = web3d_residues
                try:
                    web3d_resn_1 = aa_3_1[web3d_resn]
                except KeyError:
                    web3d_resn_1 = "."
                # Write http://weblogo.threeplusone.com/ file
                FASTA_text = "> %s %s %s %s %s, dist=%s resi=%s resn=%s %s" %(target_sel, group_id, resi, resn_1, resn_3, dist, web3d_resi, web3d_resn_1, web3d_resn)
                weblogo = "."*(var_dist)
                weblogo = weblogo[:i-1] + web3d_resn_1 + weblogo[i:]
                # Write
                wfileweblogo.write(FASTA_text + "\n")
                wfileweblogo.write(weblogo + "\n")

        # Store info
        slist.append([target_sel, group_id, resn_resi, resn_1, resi, MSA, SASA, nr_atoms_around, chain, chain_sec])

    # Group selections
    cmd.group(group, "%s_%s_*" % (target_sel, group_id))
    cmd.select("%s_sel"%group, "%s_%s_*" % (target_sel, group_id))
    # Group around
    cmd.group(group_chain, "%s_%s_%s*" % ("chain",target_sel, group_id) )
    cmd.group(group_atom, "%s_%s_%s_*" % (atom, target_sel, group_id))
    cmd.group(group_atom, "%s_around_%s_%s_*" % (atom, target_sel, group_id))
    cmd.group(group_3dweb, "%s_3dweb_sel_%s_%s_*" % (atom, target_sel, group_id))

    # Write output
    wfile = open("resi_stats_%s_%s.csv" % (target_sel, group_id), 'w')
    wfile.write("target_sel;group_id;resn_resi;resn;resi;MSA;SASA;nr_atoms_around;chain;chain_sec"+"\n")
    for i in slist:
        wfile.write("%s;%s;%s;%s;%i;%3.0f;%3.0f;%i;%s;%s" % (i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]) + "\n")
    wfile.close()
    wfileweblogo.close()

    # Back to before
    cmd.set("dot_solvent", ini_setting)
    cmd.set('dot_density', ini_setting2)

cmd.extend("match_peptides", match_peptides)
