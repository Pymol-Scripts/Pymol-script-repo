#############################################################################
#
# Author: Michel F. SANNER, Garrett MORRIS
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /opt/cvs/python/packages/share1.5/MolKit/PDBdict.py,v 1.2 2001/05/23 18:26:29 rhuey Exp $
#
# $Id: PDBdict.py,v 1.2 2001/05/23 18:26:29 rhuey Exp $
#

PDBformat = {}
PDBformat["HEADER"] = "HEADER    %.40s%.9s   %.4s" # classification, depDate, idCode
PDBformat["OBSLTE"] = "OBSLTE%.2s  %.9s %.4s %.4s      %.4s %.4s %.4s %.4s %.4s %.4s %.4s " # continuation, repDate, idCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode, rIdCode
PDBformat["TITLE "] = "TITLE   %.2s%60s" # continuation, title
PDBformat["CAVEAT"] = "CAVEAT%.2s  %.4s %s    " # continuation, idCode, comment
PDBformat["COMPND"] = "COMPND  %.2s%.60s" # continuation, compound
PDBformat["SOURCE"] = "SOURCE  %.2s%.60s" # continuation, srcName
PDBformat["KEYWDS"] = "KEYWDS  %.2s%.60s" # continuation, keywds
PDBformat["EXPDTA"] = "EXPDTA  %.2s%.60s" # continuation, technique
PDBformat["REVDAT"] = "REVDAT %3d%.2s %.9s %.5s   %1d       %.6s %.6s %.6s %.6s" # modNum, continuation, modDate, modId, modType, record, record, record, record
PDBformat["SPRSDE"] = "SPRSDE%.2s  %.9s %.4s %.4s      %.4s %.4s %.4s %.4s %.4s %.4s %.4s " # continuation, sprsdeDate, idCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode, sIdCode
PDBformat["JRNL  "] = "JRNL  %s      " # text
#PDBformat["JRNL"] = "JRNL  %.4s      %.2s%.51s " # "AUTH", continuation, authorList
#PDBformat["JRNL"] = "JRNL  %.4s      %.2s%s " # "TITL", continuation, title
#PDBformat["JRNL"] = "JRNL  %.4s      %.2s%.51s " # "EDIT", continuation, editorList
#PDBformat["JRNL"] = "JRNL  %.3s      %.15s   " # "REF", "TO BE PUBLISHED"
#PDBformat["JRNL"] = "JRNL  %.3s      %.2s%s %.2s  %s%s %4d " # "REF", continuation, pubName, "V.", volume, page, year
#PDBformat["JRNL"] = "JRNL  %.4s      %.2s%s " # "PUBL", continuation, pub
#PDBformat["JRNL"] = "JRNL  %.4s      %.4s                                                  " # "REFN", "0353"
#PDBformat["JRNL"] = "JRNL  %.4s      %.4s   %.6s %.2s  %.4s %s %.4s " # "REFN", "ASTM", astm, country, "ISBN", isbn, coden
PDBformat["REMARK"] = "REMARK %3d %s" # remarkNum, empty
#PDBformat["REMARK"] = "REMARK%.1s   %.9s %49d " # "1", "REFERENCE", refNum
#PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%.51s " # "1", "AUTH", continuation, authorList
#PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%s " # "1", "TITL", continuation, title
#PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%s " # "1", "EDIT", continuation, editorList
#PDBformat["REMARK"] = "REMARK%.1s   %.3s  %.15s   " # "1", "REF", "TO BE PUBLISHED"
#PDBformat["REMARK"] = "REMARK%.1s   %.3s  %.2s%s %.2s  %s%s %4d " # "1", "REF", continuation, pubName, "V.", volume, page, year
#PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.2s%s " # "1", "PUBL", continuation, pub
#PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.4s                                                  " # "1", "REFN", "0353"
#PDBformat["REMARK"] = "REMARK%.1s   %.4s  %.4s   %s %s  %.4s %s %.4s  " # "1", "REFN", "ASTM", astm, country, "ISBN", isbn, coden
#PDBformat["REMARK"] = "REMARK%.1s   %.11s %5.2f%.10s " # "2", "RESOLUTION.", resolution, "ANGSTROMS."
#PDBformat["REMARK"] = "REMARK%.1s   %.28s %s  " # "2", "RESOLUTION., comment
#PDBformat["REMARK"] = "REMARK%.1s   %.11s %s " # "2", "RESOLUTION.", comment
PDBformat["DBREF "] = "DBREF %.4s %.1s %4d %.1s%4d %.1s%s %s %s %5d %.1s%5d %.1s" # idCode, chainID, seqBegin, insertBegin, seqEnd, insertEnd, database, dbAccession, dbIdCode, dbseqBegin, idbnsBeg, dbseqEnd, dbinsEnd
# FIXME
#PDBformat["SEQADV"] = "SEQADV%.4s %.3s %.1s %4d %.1s%s %s %.3s %5d %s " # idCode, resName, chainID, seqNum, iCode, database, dbIdCode, dbRes, dbSeq, conflict
PDBformat["SEQRES"] = "SEQRES%2d  %.1s %4d %.3s  %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s %.3s " # serNum, chainID, numRes, resName, resName, resName, resName, resName, resName, resName, resName, resName, resName, resName, resName, resName
PDBformat["MODRES"] = "MODRES%.4s %.3s %.1s %4d %.1s%.3s %s  " # idCode, resName, chainID, seqNum, iCode, stdRes, comment
PDBformat["HET   "] = "HET   %.3s %.1s  %4d%.1s%5d  %s     " # hetID, ChainID, seqNum, iCode, numHetAtoms, text
PDBformat["HETNAM"] = "HETNAM%.2s  %.3s %s " # continuation, hetID, text
PDBformat["HETSYN"] = "HETSYN%.2s  %.3s %.55s " # continuation, hetID, hetSynonyms
PDBformat["FORMUL"] = "FORMUL%2d  %.3s  %2d %.1s%s" # compNum, hetID, continuation, asterisk, text
PDBformat["HELIX "] = "HELIX %3d %.3s %.3s %.1s %4d %.1s%.3s %.1s %4d %.1s%2d%s%5d " # serNum, helixID, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum, endICode, helixClass, comment, length
PDBformat["SHEET "] = "SHEET %3d %.3s %2d%.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%2d%.4s %.3s%.1s %4d%.1s%.4s %.3s%.1s %4d%.1s" # strand, sheetID, numStrands, initResName, initChainID, initSeqNum, initICode, endResName, endChainID, endSeqNum, endICode, sense, curAtom, curResName, curChainId, curResSeq, curICode, prevAtom, prevResName, prevChainId, prevResSeq, prevICode
PDBformat["TURN  "] = "TURN  %3d %.3s %.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%s    " # seq, turnId, initResName, initChainId, initSeqNum, initICode, endResName, endChainId, endSeqNum, endICode, comment
PDBformat["SSBOND"] = "SSBOND%3d %.3s %.1s %4d %.1s%.3s   %.1s %4d %.1s%.6s                       %.6s " # serNum, "CYS", chainID1, seqNum1, icode1, "CYS", chainID2, seqNum2, icode2, sym1, sym2
PDBformat["LINK  "] = "LINK  %.4s      %.1s%.3s%.1s %4d%.1s%.4s               %.1s%.3s%.1s %4d%.1s%.6s  %.6s " # name1, altLoc1, resName1, chainID1, resSeq1, iCode1, name2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1, sym2
PDBformat["HYDBND"] = "HYDBND%.4s      %.1s%.3s%.1s %5d%.1s%.4s %.1s%.1s %5d%.1s%.4s %.1s%.3s%.1s %5d%.1s%.6s%.6s " # name1, altLoc1, resName1, Chain1, resSeq1, ICode1, nameH, altLocH, ChainH, resSeqH, iCodeH, name2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1, sym2
PDBformat["SLTBRG"] = "SLTBRG%.4s      %.1s%.3s%.1s %4d%.1s%.4s               %.1s%.3s%.1s %4d%.1s%.6s  %.6s " # atom1, altLoc1, resName1, chainID1, resSeq1, iCode1, atom2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1, sym2
PDBformat["CISPEP"] = "CISPEP%3d %.3s %.1s %4d %.1s%.3s   %.1s %4d %.1s%3d       %6.2f       " # serNum, pep1, chainID1, seqNum1, icode1, pep2, chainID2, seqNum2, icode2, modNum, measure
PDBformat["SITE  "] = "SITE  %3d %.3s %2d %.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%.3s %.1s %4d%.1s%.3s %.1s %4d%.1s" # seqNum, siteID, numRes, resName1, chainID1, seq1, iCode1, resName2, chainID2, seq2, iCode2, resName3, chainID3, seq3, iCode3, resName4, chainID4, seq4, iCode4
PDBformat["CRYST1"] = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s %4d" # a, b, c, alpha, beta, gamma, sGroup, z
PDBformat["ORIGXn"] = "ORIGXn%10.6f    %10.6f%10.6f%10.5f     " # o[n][1], o[n][2], o[n][3], t[n]
PDBformat["SCALEn"] = "SCALEn%10.6f    %10.6f%10.6f%10.5f     " # s[n][1], s[n][2], s[n][3], u[n]
PDBformat["MTRIXn"] = "MTRIXn%3d %10.6f%10.6f%10.6f%10.5f     %1d    " # serial, m[n][1], m[n][2], m[n][3], v[n], iGiven
PDBformat["TVECT "] = "TVECT %3d %10.5f%10.5f%10.5f%s" # serial, t[1], t[2], t[3], text
PDBformat["MODEL "] = "MODEL %4d    " # serial
PDBformat["ATOM  "] = "ATOM  %5d%.4s %.1s%.3s%.1s %4d%.1s%8.3f   %8.3f%8.3f%6.2f%6.2f%.4s      %.2s%.2s" # serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, segID, element, charge
PDBformat["SIGATM"] = "SIGATM%5d%.4s %.1s%.3s%.1s %4d%.1s%8.3f   %8.3f%8.3f%6.2f%6.2f%.4s      %.2s%.2s" # serial, name, altLoc, resName, chainID, resSeq, iCode, sigX, sigY, sigZ, sigOcc, sigTemp, segID, element, charge
PDBformat["ANISOU"] = "ANISOU%5d%.4s %.1s%.3s%.1s %4d%.1s%7d %7d%7d%7d%7d%7d%.4s  %.2s%.2s" # serial, name, altLoc, resName, chainID, resSeq, iCode, u[0][0], u[1][1], u[2][2], u[0][1], u[0][2], u[1][2], segID, element, charge
PDBformat["SIGUIJ"] = "SIGUIJ%5d%.4s %.1s%.3s%.1s %4d%.1s%7d %7d%7d%7d%7d%7d%.4s  %.2s%.2s" # serial, name, altLoc, resName, chainID, resSeq, iCode, sig[1][1], sig[2][2], sig[3][3], sig[1][2], sig[1][3], sig[2][3], segID, element, charge
PDBformat["TER   "] = "TER   %5d%.3s      %.1s %4d%.1s" # serial, resName, chainID, resSeq, iCode
PDBformat["HETATM"] = "HETATM%5d%.4s %.1s%.3s%.1s %4d%.1s%8.3f   %8.3f%8.3f%6.2f%6.2f%.4s      %.2s%.2s" # serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, segID, element, charge
PDBformat["ENDMDL"] = "ENDMDL" # 
PDBformat["CONECT"] = "CONECT%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d" # serial, serial, serial, serial, serial, serial, serial, serial, serial, serial, serial
PDBformat["MASTER"] = "MASTER%5d    %5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d" # numRemark, 0, numHet, numHelix, numSheet, numTurn, numSite, numXform, numCoord, numTer, numConect, numSeq


PDBFormatConstr = {}
PDBFormatConstr["SSBOND"] = [3, None, None, 4, None, None, None, 4, None, None, None]
PDBFormatConstr["REVDAT"] = [3, None, None, None, 1, None, None, None, None]
PDBFormatConstr["MODEL "] = [4]
PDBFormatConstr["SCALEn"] = [10, 10, 10, 10]
PDBFormatConstr["COMPND"] = [None, None]
PDBFormatConstr["SPRSDE"] = [None, None, None, None, None, None, None, None, None, None, None]
PDBFormatConstr["CAVEAT"] = [None, None, None]
PDBFormatConstr["HETNAM"] = [None, None, None]
PDBFormatConstr["TVECT "] = [3, 10, 10, 10, None]
PDBFormatConstr["FORMUL"] = [2, None, 2, None, None]
PDBFormatConstr["SOURCE"] = [None, None]
PDBFormatConstr["MODRES"] = [None, None, None, 4, None, None, None]
PDBFormatConstr["HETATM"] = [5, None, None, None, None, 4, None, 8, 8, 8, 6, 6, None, None, None]
PDBFormatConstr["ENDMDL"] = []
PDBFormatConstr["ATOM  "] = [5, None, None, None, None, 4, None, 8, 8, 8, 6, 6, None, None, None]
PDBFormatConstr["CISPEP"] = [3, None, None, 4, None, None, None, 4, None, 3, 6]
PDBFormatConstr["CRYST1"] = [9, 9, 9, 7, 7, 7, None, 4]
PDBFormatConstr["SLTBRG"] = [None, None, None, None, 4, None, None, None, None, None, 4, None, None, None]
PDBFormatConstr["EXPDTA"] = [None, None]
PDBFormatConstr["SIGUIJ"] = [5, None, None, None, None, 4, None, 7, 7, 7, 7, 7, 7, None, None, None]
PDBFormatConstr["HET   "] = [None, None, 4, None, 5, None]
PDBFormatConstr["SHEET "] = [3, None, 2, None, None, 4, None, None, None, 4, None, 2, None, None, None, 4, None, None, None, None, 4, None]
PDBFormatConstr["HYDBND"] = [None, None, None, None, 5, None, None, None, None, 5, None, None, None, None, None, 5, None, None, None]
PDBFormatConstr["REMARK"] = [3, None]
PDBFormatConstr["TITLE "] = [None, 60]
PDBFormatConstr["ANISOU"] = [5, None, None, None, None, 4, None, 7, 7, 7, 7, 7, 7, None, None, None]
PDBFormatConstr["SIGATM"] = [5, None, None, None, None, 4, None, 8, 8, 8, 6, 6, None, None, None]
PDBFormatConstr["TURN  "] = [3, None, None, None, 4, None, None, None, 4, None, None]
PDBFormatConstr["SITE  "] = [3, None, 2, None, None, 4, None, None, None, 4, None, None, None, 4, None, None, None, 4, None]
PDBFormatConstr["HELIX "] = [3, None, None, None, 4, None, None, None, 4, None, 2, None, 5]
PDBFormatConstr["TER   "] = [5, None, None, 4, None]
PDBFormatConstr["HEADER"] = [None, None, None]
PDBFormatConstr["MASTER"] = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
PDBFormatConstr["CONECT"] = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
PDBFormatConstr["HETSYN"] = [None, None, None]
PDBFormatConstr["DBREF "] = [None, None, 4, None, 4, None, None, None, None, 5, None, 5, None]
PDBFormatConstr["MTRIXn"] = [3, 10, 10, 10, 10, 1]
PDBFormatConstr["ORIGXn"] = [10, 10, 10, 10]
PDBFormatConstr["JRNL  "] = [None]
PDBFormatConstr["SEQRES"] = [2, None, 4, None, None, None, None, None, None, None, None, None, None, None, None, None]
PDBFormatConstr["LINK  "] = [None, None, None, None, 4, None, None, None, None, None, 4, None, None, None]
PDBFormatConstr["OBSLTE"] = [None, None, None, None, None, None, None, None, None, None, None]
PDBFormatConstr["KEYWDS"] = [None, None]
















