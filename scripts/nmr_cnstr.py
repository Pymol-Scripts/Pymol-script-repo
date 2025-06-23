'''
See more here: http://www.pymolwiki.org/index.php/nmr_cnstr
##############################################################################################################
# Pymol Script: For visualizing the NMR constrains (DYANA & CNS format), on top of the calculated structure. #
#               Author: Evangelos Papadopoulos.                                                              #
#  previous affiliation: Dept. of Biochemistry and Biophysics,                                               #
#                       Arrhenius Laboratories,                                                              #
#                       Stockholm University                                                                 #
#                       SE-106 91 Stockholm, Sweden                                                          #
#                email:evangelos.papadopoulos@gmail.com                                                      #
#                NOTES: This is a preliminary version.                                                       #
#                                                                                                            #
#     Reference: please if you find this script usefull add the following reference:                   #
#     NMR Solution Structure of the Peptide Fragment 1-30, Derived from Unprocessed Mouse Doppel             #
#     Protein, in DHPC Micelles. Biochemistry. 2006 Jan 10;45(1):159-166. PMID: 16388591                     #
#                                                                                                            #
##############################################################################################################
'''

from __future__ import print_function

def upl(fname):

    f = open(fname, 'r')
    i = 1
    upl = f.readline()
#
    while upl != '':

        print(upl, i)
        cns = upl.split()
        cmd.dist('upl' + str(i), 'i. ' + cns[0] + ' & n. ' + cns[2], 'i. ' + cns[3] + ' & n. ' + cns[5])
        upl = f.readline()
        i += 1
#
    f.close()
    cmd.hide('labels')
    cmd.set('dash_gap', 0.05)
    cmd.do("orient")
    cmd.set('movie_delay', 1500)


def cns(fname):

    f = open(fname, 'r')
    i = 1
    upl = f.readline()
    print(upl, i)
    while upl != '':
        if upl == '\n':
            upl = f.readline()
            continue
        cns = upl.split()
        print(cns, i)
        if cns[0] == 'assign':
            print('CNS')
            if cns[5] == 'HB*':
                print('CNS***')
            cmd.dist('upl' + str(i), 'i. ' + cns[2] + ' & n. ' + cns[5], 'i. ' + cns[7] + ' & n. ' + cns[10])
        i += 1
        upl = f.readline()
        print('*' + upl + '*', i)

    f.close()
    cmd.set('dash_gap', 0.05)
    cmd.hide('labels')
    cmd.do("orient")
    cmd.set('movie_delay', 1500)
