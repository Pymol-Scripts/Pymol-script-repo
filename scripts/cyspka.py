from __future__ import print_function
from pymol import cmd
import os
import sys
import math
from time import localtime, strftime

# Thx for inspiration from Simple scriptin PymMOl http://www.pymolwiki.org/index.php/Simple_Scripting
# Made by Ma.Sc student. Troels Linnet, 2011-08. troels.linnet@bbz.uni-leipzig.de
# Based solely on the work by:
# Maik H. Jacob, Dan Amir, Vladimir Ratner, Eugene Gussakowsky, and Elisha Haas
# Predicting Reactivities of Protein Surface Cysteines as Part of a Strategy for Selective Multiple Labeling. (Biochemistry 2005, 44, 13664-13672)

# Example of pymol script: Directory "predict_reactivity" has script file cyspka.py and cysteine residue pdb file: cys.pdb
# import cyspka
# fetch 4AKE, async=0
# create 4AKE-A, /4AKE//A and not resn HOH
# delete 4AKE
# hide everything
# show cartoon, 4AKE-A
# cyspka 4AKE-A, A, 18


def cyspka(molecule, chain, residue, SeeProgress='yes', pH=7.2, MoveSGatom='no', SGatom=str((0, 0, 0))):
    # If SeeProgress='yes', computation time will take 10-20% extra, but nice to follow.
    cmd.refresh()
    RotationRange = 360
    RotationDegree = 1
    # For error checking, the energies can be printed out
    printMC = 'no'
    printSC = 'no'
    # Parameters
    DieElecSpheDist = 7.0
    DieElecWaterDist = 1.4
    DieElecWater = 78.5
    DieElecCore = 4.0
    BornPenaltyB = 1.0
    AvogadroR = 8.31446216
    Temp = 298
    DeltapKMCSC = 0
    pK1 = 9.25
    pK2 = 8.0
    NotPopuDist = 2.4
    PopEnergyPenalty = 10000000
    # Side chain discrete charges
    DieElecSC = 40.0
    SCchargeASP = -1
    SCchargeGLU = -1
    SCchargeOXT = -1
    SCchargeARG = +1
    SCchargeHIS = +1
    SCchargeLYS = +1
    SCchargeMET1 = +1
    # Main chain partial charges
    NrMainchainNeighBours = 5
    DieElecMC = 22.0
    MCchargeC = +0.55
    MCchargeO = -0.55
    MCchargeN = -0.35
    MCchargeH = +0.35
    MCchargeProCA = +0.1
    MCchargeProCD = +0.1
    MCchargeProN = -0.2

    # Loading an Cys residue, give it a logic name, and aligning it. The oxygen atom can not be aligned in many cases, and are skipped.
    # We use only this molecule, to find the initial position of the SG atom, and to rotate the SG atom around the CA-CB bond. The molecule atom positions are not used for electric potential calculatons.
    Cysmolecule = str(molecule) + str(residue) + "Cys"
    cmd.fragment("cys")
    cmd.set_name('cys', Cysmolecule)
    # We use pair_fir, since align and super gets unstable with so few atoms
    pairfitCys(Cysmolecule, molecule, chain, residue)
    # Give nice representations quickly
    cmd.show("sticks", Cysmolecule)
    cmd.select(str(molecule) + str(residue) + "Res", "/" + molecule + "//" + chain + "/" + residue)
    print("/" + molecule + "//" + chain + "/" + residue)
    cmd.show("sticks", str(molecule) + str(residue) + "Res")
    cmd.disable(str(molecule) + str(residue) + "Res")
    # Find out what is the residuename we are investigating for
    Respdbstr = cmd.get_pdbstr(str(molecule) + str(residue) + "Res")
    Ressplit = Respdbstr.split()
    residueName = Ressplit[3]

    print("")
    print("# Hello, PyMOLers. It should take around 1 minute per residue.")
    print("# molecule: %s , chain: %s, residue: %s %s, pH: %s " % (molecule, chain, residueName, residue, pH))

    # Determine the range of neighbour residues possible.
    Maxresidues = cmd.count_atoms("/" + molecule + "//" + chain + " and name CA")
    for i in range(NrMainchainNeighBours + 1):
        if int(residue) - i >= 1:
            Minresidue = int(residue) - i
        else:
            break
    for i in range(NrMainchainNeighBours + 1):
        if int(residue) + i <= Maxresidues:
            Maxresidue = int(residue) + i
        else:
            break

    # Get the position and the vector for the CA->CB bond.
    dihedN = "/" + Cysmolecule + "//" + "/" + "/N"
    dihedCA = "/" + Cysmolecule + "//" + "/" + "/CA"
    dihedCB = "/" + Cysmolecule + "//" + "/" + "/CB"
    dihedSG = "/" + Cysmolecule + "//" + "/" + "/SG"
    dihedralPosCA = cmd.get_atom_coords(dihedCA)
    dihedralPosSG = cmd.get_atom_coords(dihedSG)
    dihedralVector = AtomVector(dihedCA, dihedCB)

    # To compare with article, we can move the SGatom to a starting position. The rotation is still determined around the CA-CB bond.
    if MoveSGatom == 'yes':
        SGatom = [float(SGatom[1:-1].split(",")[0]), float(SGatom[1:-1].split(",")[1]), float(SGatom[1:-1].split(",")[2])]
        Translate = [(SGatom[0] - dihedralPosSG[0]), (SGatom[1] - dihedralPosSG[1]), (SGatom[2] - dihedralPosSG[2])]
        cmd.translate(Translate, dihedSG, state=0, camera=0)
        dihedralPosSG = cmd.get_atom_coords(dihedSG)
    # Create a pymol molecule, that in the end will hold and show all SG atoms. Gives the representation of the rotameric states.
    SGName = str(molecule) + str(residue) + "SG"
    cmd.create(SGName, "None")
    # Create a pymol molecule, that in the end will hold and show all Amide protons. Gives a nice representation, and easy to delete.
    AmideName = str(molecule) + str(residue) + "NH"
    cmd.create(AmideName, "None")
    # Check if there are any nearby SG atoms, which could make a SG-SG dimer formation. The
    breakDimer = "no"
    breakDimer = CheckDimer(dihedSG, molecule, chain, residue)
    # Create a list for appending the calculated energies.
    ListofEnergies = []
    ListofRotamerDiscarded = []
    # print "Angle before rotation", cmd.get_dihedral(dihedN,dihedCA,dihedCB,dihedSG)
    # Enter into the loop of rotameric states
    for i in range(int(math.floor(RotationRange / RotationDegree))):
        Angle = i * RotationDegree
        # Create pymol molecule/SG atom for which we will calculate for.
        SGNameAngle = str(residue) + "SG" + str(Angle)
        cmd.create(SGNameAngle, dihedSG)
        # Calculate new coordinates for rotation around CA->CB bond. Then translate the created SG atom.
        SGNewPos = fRotateAroundLine(dihedralPosSG, dihedralPosCA, dihedralVector, Angle)
        Translate = [(SGNewPos[0] - dihedralPosSG[0]), (SGNewPos[1] - dihedralPosSG[1]), (SGNewPos[2] - dihedralPosSG[2])]
        cmd.translate(Translate, SGNameAngle, state=0, camera=0)
        # If one wants to "see it happen" while its making the states. But it will take extra computation time.
        if SeeProgress == 'yes':
            cmd.refresh()
        # Calculate number of neighbours within 2.4 Angstrom. Amide hydrogens are not considered, and are actually not build yet.
        nameselect = "(((/" + molecule + "//" + chain + " and not /" + molecule + "//" + chain + "/" + residue + ") or /" + molecule + "//" + chain + "/" + residue + "/N+CA+C+O)  within " + str(NotPopuDist) + " of /" + SGNameAngle + "//" + "/" + "/SG) and not resn HOH"
        # print nameselect
        cmd.select("NotPop", nameselect)
        NotPopNr = cmd.count_atoms("NotPop")
        # print Angle, NotPopNr, cmd.get_dihedral(dihedN,dihedCA,dihedCB,SGNameAngle)
        # If no neighbours, then proceed calculating
        if NotPopNr == 0:
            SumAllWMC = 0.0
            # Now calculate the electric potential due to the side chains.
            SumWSC = fSumWSC(molecule, SGNameAngle, chain, residue, DieElecSC, SCchargeASP, SCchargeGLU, SCchargeOXT, SCchargeARG, SCchargeHIS, SCchargeLYS, SCchargeMET1, printSC)
            # Now we calculate for the flanking 5 peptide groups on each side of the Cysteine CA atom.
            # For the first residue, only calculate for the tailing C,O atom in the peptide bond. No test for Proline.
            SumWMCFirst = fSumWMCFirst(molecule, SGNameAngle, chain, residue, Minresidue, DieElecMC, MCchargeC, MCchargeO, printMC)
            # For the residue itself, we dont test for PRO, since it should be a Cysteine.
            SumWMCresidue = fSumWMCresidue(molecule, SGNameAngle, chain, residue, int(residue), DieElecMC, MCchargeC, MCchargeO, MCchargeN, MCchargeH, AmideName, printMC)
            # For the last residue, we test for Proline. We only calculate for the N,H atom, or if Proline, N,CA and CD atom.
            SumWMCLast = fSumWMCLast(molecule, SGNameAngle, chain, residue, Maxresidue, DieElecMC, MCchargeN, MCchargeH, MCchargeProCA, MCchargeProCD, MCchargeProN, AmideName, printMC)
            # Then loop over rest of the residues in the chain.
            for j in (list(range(Minresidue + 1, int(residue))) + list(range(int(residue) + 1, Maxresidue))):
                MCNeighbour = j
                # print "Looking at neighbour", j
                SumWMC = fSumWMC(molecule, SGNameAngle, chain, residue, MCNeighbour, DieElecMC, MCchargeC, MCchargeO, MCchargeN, MCchargeH, MCchargeProCA, MCchargeProCD, MCchargeProN, AmideName, printMC)
                SumAllWMC = SumAllWMC + SumWMC
                # print "Rotation: %s Neighbour: %s " % (Angle, j)
            # Since the SG atom is negative, we multiply with -1.
            SumMCSC = -1 * (SumWSC + SumWMCFirst + SumWMCresidue + SumWMCLast + SumAllWMC)
            # Makes the neighbour count. Everything in 'molecule" within 7 ang of aligned SG atom. Not counting 'residue'. Adding 5 for 'residue' N,CA,C,O,CB
            ListNeighbourCount = fNeighbourCount(molecule, SGNameAngle, chain, residue, DieElecSpheDist)
            # Calculate the weighted electric potential and alter the b factor for coloring. Then add the rotated SG into bucket of SG atoms.
            SG_MCSC_Weight = fBoltzSingleState(SumMCSC, AvogadroR, Temp) * SumMCSC
            cmd.alter(SGNameAngle, 'b="%s"' % SG_MCSC_Weight)
            cmd.alter(SGNameAngle, 'name="S%s"' % Angle)
            cmd.create(SGName, SGName + " + " + SGNameAngle)
            # Then save the calculated values
            ListofEnergies.append([Angle, SumMCSC, ListNeighbourCount, NotPopNr, SG_MCSC_Weight, cmd.get_atom_coords(SGNameAngle)])
            cmd.delete(SGNameAngle)
        else:
            SumMCSCPenalty = PopEnergyPenalty
            ListNeighbourCount = fNeighbourCount(molecule, SGNameAngle, chain, residue, DieElecSpheDist)
            ListofRotamerDiscarded.append([Angle, SumMCSCPenalty, ListNeighbourCount, NotPopNr, 0, cmd.get_atom_coords(SGNameAngle)])
            cmd.delete(SGNameAngle)
    # Now show all the SG atoms as the available rotameric states.
    cmd.show("nb_spheres", SGName)
    cmd.delete("NotPop")
    cmd.spectrum("b", selection=SGName)
    AvailRotStates = len(ListofEnergies)
    # print "Available Rotational States: ", AvailRotStates

    # Do the calculations according to eq 5.
    # Find the partition function
    BoltzPartition = 0.0
    for i in range(len(ListofEnergies)):
        Boltz = fBoltzSingleState(ListofEnergies[i][1], AvogadroR, Temp)
        BoltzPartition = BoltzPartition + Boltz
    # Find the summed function
    BoltzSumNi = 0.0
    for i in range(len(ListofEnergies)):
        BoltzNi = fBoltzSingleState(ListofEnergies[i][1], AvogadroR, Temp) * ListofEnergies[i][1]
        BoltzSumNi = BoltzSumNi + BoltzNi

    # Check if there was any possible rotamers

    nostates = "no"
    if len(ListofEnergies) == 0:
        print("####################################################")
        print("########### WARNING: No states available ###########")
        print("########### Did you mutate a Glycine?    ###########")
        print("####################################################")
        BoltzSumNi = 0
        BoltzPartition = 0
        BoltzMCSC = 0
        DeltapKMCSC = 99
        NeighbourCount = 0
        nostates = "yes"
    else:
        # Final calculation
        BoltzMCSC = (BoltzSumNi) / (BoltzPartition)
        DeltapKMCSC = fDeltapK(BoltzMCSC, AvogadroR, Temp)

    # Find average number of neighbours
    NCSum = 0.0
    NCWeightedSum = 0.0
    for i in range(len(ListofEnergies)):
        NCi = ListofEnergies[i][2]
        NCSum = NCSum + NCi
        NCWeightedi = fBoltzSingleState(ListofEnergies[i][1], AvogadroR, Temp) * ListofEnergies[i][2] / BoltzPartition
        NCWeightedSum = NCWeightedSum + NCWeightedi
    # print "Weighted neighbour", int(round(NCWeightedSum))
    #NeighbourCount = int(round(NCSum/len(ListofEnergies)))
        NeighbourCount = round(NCWeightedSum, 1)
    # If we found dimers
    if breakDimer == "yes":
        print("####################################################")
        print("########### WARNING: Dimer formation?    ###########")
        print("####################################################")
        BoltzSumNi = 0
        BoltzPartition = 0
        BoltzMCSC = 0
        DeltapKMCSC = 99
        NeighbourCount = 0

    # Calculate the BornPenalty based on the neighbour count. It's a wrapper script for equation 13, 12, 11.
    EnergyBornPenalty = fEnergyBornPenalty(DieElecSpheDist, DieElecWaterDist, NeighbourCount, DieElecWater, DieElecCore, BornPenaltyB)
    DeltapKB = fDeltapK(EnergyBornPenalty, AvogadroR, Temp)

    # Do the calculations according to eq 3 and 9.
    pKm1 = fpKm1(DeltapKMCSC, pK1)
    pKm2 = fpKm2(DeltapKMCSC, DeltapKB, pK2)
    FracCysm1 = fFracCys(pKm1, pH)
    FracCysm2 = fFracCys(pKm2, pH)

    # Lets make a result file, and write out the angle, the SumMCSC, and the number of neighbours for this state.
    Currentdir = os.getcwd()
    Newdir = os.path.join(os.getcwd(), "Results")
    if not os.path.exists(Newdir):
        os.makedirs(Newdir)
    filename = os.path.join(".", "Results", "Result_" + molecule + "_" + chain + "_" + residue + ".txt")
    filenamelog = os.path.join(".", "Results", "Result_log.log")
    logfile = open(filenamelog, "a")
    outfile = open(filename, "w")
    timeforlog = strftime("%Y %b %d %a %H:%M:%S", localtime())
    logfile.write("# " + timeforlog + "\n")
    logfile.write("# molecule: %s , chain: %s, residue: %s %s, pH: %s " % (molecule, chain, residueName, residue, pH) + "\n")
    logfile.write("# BoltzSumNi:  BoltzPartition:  BoltzMCSC" + "\n")
    logfile.write(("# %.4f  %.4f  %.4f" + '\n') % (BoltzSumNi, BoltzPartition, BoltzMCSC))
    logfile.write("#    Res  NC    States  pKmcsc  pK1   pKB     pK2  pKm1     pKm2    f(C-)m1   f(C-)m2" + "\n")
    logfile.write(("; %s %s   %s  %s     %.4f  %s  %.4f  %s  %.4f  %.4f  %.6f  %.6f" + '\n') % (residueName, residue, NeighbourCount, AvailRotStates, DeltapKMCSC, pK1, DeltapKB, pK2, pKm1, pKm2, FracCysm1, FracCysm2))
    if nostates == "yes":
        logfile.write("##### ERROR; No states available ###" + "\n")
    if breakDimer == "yes":
        logfile.write("##### ERROR; Dimer formation ###" + "\n")
    logfile.write('\n')
    outfile.write("# molecule: %s , chain: %s, residue: %s %s, pH: %s " % (molecule, chain, residueName, residue, pH) + "\n")
    outfile.write("# BoltzSumNi:  BoltzPartition:  BoltzMCSC" + "\n")
    outfile.write(("# %.4f  %.4f  %.4f" + '\n') % (BoltzSumNi, BoltzPartition, BoltzMCSC))
    outfile.write("#    Res  NC    States  pKmcsc  pK1   pKB     pK2  pKm1     pKm2    f(C-)m1   f(C-)m2" + "\n")
    outfile.write(("; %s %s   %s  %s     %.4f  %s  %.4f  %s  %.4f  %.4f  %.6f  %.6f" + '\n') % (residueName, residue, NeighbourCount, AvailRotStates, DeltapKMCSC, pK1, DeltapKB, pK2, pKm1, pKm2, FracCysm1, FracCysm2))
    if nostates == "yes":
        outfile.write("##### ERROR; No states available ###" + "\n")
    if breakDimer == "yes":
        outfile.write("##### ERROR; Dimer formation ###" + "\n")
    outfile.write('\n')
    outfile.write("#Ang  SumMCSC   NC rNC MCSC_Weight       SG[X,Y,Z]" + "\n")
    for i in range(len(ListofEnergies)):
        outfile.write("%4.1d %10.3f %2.1d %1.1d %10.3f [%8.3f, %8.3f, %8.3f]" % (ListofEnergies[i][0], ListofEnergies[i][1], ListofEnergies[i][2], ListofEnergies[i][3], ListofEnergies[i][4], ListofEnergies[i][5][0], ListofEnergies[i][5][1], ListofEnergies[i][5][2]) + '\n')
    for i in range(len(ListofRotamerDiscarded)):
        outfile.write("%4.1d %10.3f %2.1d %1.1d %10.3f [%8.3f, %8.3f, %8.3f]" % (ListofRotamerDiscarded[i][0], ListofRotamerDiscarded[i][1], ListofRotamerDiscarded[i][2], ListofRotamerDiscarded[i][3], ListofRotamerDiscarded[i][4], ListofRotamerDiscarded[i][5][0], ListofRotamerDiscarded[i][5][1], ListofRotamerDiscarded[i][5][2]) + '\n')
    outfile.close()

    # Now, we are done. Just print out. The ; is for a grep command to select these lines in the output.
    print("# residue: %s %s. Average NeighbourCount NC= %s " % (residueName, residue, NeighbourCount))
    print("# From residue %s to residue %s" % (Minresidue, Maxresidue))
    print("# BoltzSumNi:  BoltzPartition:  BoltzMCSC")
    print("# %.4f  %.4f  %.4f" % (BoltzSumNi, BoltzPartition, BoltzMCSC))
    print("# Result written in file: %s" % (filename))
    print("#    Res  NC    States  pKmcsc  pK1   pKB     pK2  pKm1     pKm2    f(C-)m1   f(C-)m2")
    print("; %s %s   %s  %s     %.4f  %s  %.4f  %s  %.4f  %.4f  %.6f  %.6f" % (residueName, residue, NeighbourCount, AvailRotStates, DeltapKMCSC, pK1, DeltapKB, pK2, pKm1, pKm2, FracCysm1, FracCysm2))
    if nostates == "yes":
        print("##### ERROR; No states available ###")
    if breakDimer == "yes":
        print("##### ERROR; Dimer formation ###")
cmd.extend("cyspka", cyspka)


def loopcyspka(molecule, chain, residue, SeeProgress='no', pH=7.2, MoveSGatom='no', SGatom=str((0, 0, 0))):
    residue = residue.split('.')
    residueList = []
    for i in residue:
        if '-' in i:
            tmp = i.split('-')
            residueList.extend(list(range(int(tmp[0]), int(tmp[-1]) + 1)))
        if '-' not in i:
            residueList.append(int(i))
    print("Looping over residues")
    print(residueList)
    for i in residueList:
        cyspka(molecule, chain, str(i), SeeProgress, pH, MoveSGatom, SGatom)
cmd.extend("loopcyspka", loopcyspka)


def fNeighbourCount(molecule, Cysmolecule, chain, residue, DieElecSpheDist):
    nameselect = "(((/" + molecule + "//" + chain + " and not /" + molecule + "//" + chain + "/" + residue + ") or /" + molecule + "//" + chain + "/" + residue + "/N+CA+C+O)  within " + str(DieElecSpheDist) + " of /" + Cysmolecule + "//" + "/" + "/SG) and not resn HOH"
    # print nameselect
    cmd.select(residue + "NC", nameselect)
    # Adding 1 for CB
    Neighbours = cmd.count_atoms(residue + "NC") + 1
    cmd.delete(residue + "NC")
    return Neighbours


def fNeighbourWater(DieElecSpheDist, DieElecWaterDist, NeighbourCount):
    Waters = 0.74 * math.pow(DieElecSpheDist, 3) / math.pow(DieElecWaterDist, 3) - NeighbourCount
    return Waters


def fDieElecEF(NeighbourWater, DieElecWater, NeighbourCount, DieElecCore):
    DieElecEF = (NeighbourWater * DieElecWater + NeighbourCount * DieElecCore) / (NeighbourWater + NeighbourCount)
    return DieElecEF


def fBornPenalty(BornPenaltyB, DieElecEF, DieElecWater):
    BornPenalty = (1.39 * math.pow(10, 6)) / (2 * BornPenaltyB) * (1.0 / DieElecEF - 1.0 / DieElecWater)
    return BornPenalty


def fEnergyBornPenalty(DieElecSpheDist, DieElecWaterDist, NeighbourCount, DieElecWater, DieElecCore, BornPenaltyB):
    NeighbourWater = fNeighbourWater(DieElecSpheDist, DieElecWaterDist, NeighbourCount)
    DieElecEF = fDieElecEF(NeighbourWater, DieElecWater, NeighbourCount, DieElecCore)
    BornPenalty = fBornPenalty(BornPenaltyB, DieElecEF, DieElecWater)
    return BornPenalty


def fDeltapK(Energy, AvogadroR, Temp):
    DeltapK = -1 * math.log10(math.exp(-Energy / (AvogadroR * Temp)))
    return DeltapK


def fRotateAroundLine(OriPoint, ThroughLinePoint, LineVector, AngleDeg):
    # See http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/. Section 6.1
    AngleRad = math.radians(AngleDeg)
    x = OriPoint[0]
    y = OriPoint[1]
    z = OriPoint[2]
    a = ThroughLinePoint[0]
    b = ThroughLinePoint[1]
    c = ThroughLinePoint[2]
    u = LineVector[0]
    v = LineVector[1]
    w = LineVector[2]
    L = math.pow(u, 2) + math.pow(v, 2) + math.pow(w, 2)
    xPos = ((a * (math.pow(v, 2) + math.pow(w, 2)) - u * (b * v + c * w - u * x - v * y - w * z)) * (1 - math.cos(AngleRad)) + L * x * math.cos(AngleRad) + math.sqrt(L) * (-c * v + b * w - w * y + v * z) * math.sin(AngleRad)) / L
    yPos = ((b * (math.pow(u, 2) + math.pow(w, 2)) - v * (a * u + c * w - u * x - v * y - w * z)) * (1 - math.cos(AngleRad)) + L * y * math.cos(AngleRad) + math.sqrt(L) * (c * u - a * w + w * x - u * z) * math.sin(AngleRad)) / L
    zPos = ((c * (math.pow(u, 2) + math.pow(v, 2)) - w * (a * u + b * v - u * x - v * y - w * z)) * (1 - math.cos(AngleRad)) + L * z * math.cos(AngleRad) + math.sqrt(L) * (-b * u + a * v - v * x + u * y) * math.sin(AngleRad)) / L
    NewPos = [xPos, yPos, zPos]
    return NewPos


def fWSC(charge, DieElecSC, DistR):
    # print charge, DistR
    WSC = 1.39 * math.pow(10, 6) * charge / (DieElecSC * DistR)
    return WSC


def fSumWSC(molecule, SGNameAngle, chain, residue, DieElecSC, SCchargeASP, SCchargeGLU, SCchargeOXT, SCchargeARG, SCchargeHIS, SCchargeLYS, SCchargeMET1, printSC):
    SumWSC = 0.0
    SGnameselect = "/" + SGNameAngle + "//" + "/" + "/SG"
    # Sidechain ASP
    nameselect = "/" + molecule + " and resn ASP and name CG and not resi " + residue
    cmd.select("SC", nameselect)
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distASP', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeASP, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC ASP ", str(SClist[i]), " ", SCchargeASP, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distASP')
    # Sidechain GLU
    nameselect = "/" + molecule + " and resn GLU and name CD and not resi " + residue
    cmd.select("SC", nameselect)
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distGLU', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeGLU, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC GLU ", str(SClist[i]), " ", SCchargeGLU, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distGLU')
    # print "GLU", cmd.count_atoms("SC"), SumWSC
    # Sidechain OXT
    nameselect = "/" + molecule + " and byres name OXT and not resi " + residue
    cmd.select("SC", nameselect)
    cmd.select("SC", "SC and name C")
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distOXT', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeOXT, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC OXT ", str(SClist[i]), " ", SCchargeOXT, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distOXT')
    # print "OXT", cmd.count_atoms("SC"), SumWSC
    # Sidechain ARG
    nameselect = "/" + molecule + " and resn ARG and name CZ and not resi " + residue
    cmd.select("SC", nameselect)
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distARG', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeARG, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC ARG ", str(SClist[i]), " ", SCchargeARG, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distARG')
    # print "ARG", cmd.count_atoms("SC"), SumWSC
    # Sidechain HIS
    nameselect = "/" + molecule + " and resn HIS and name CD2 and not resi " + residue
    cmd.select("SC", nameselect)
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distHIS', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeHIS, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC HIS ", str(SClist[i]), " ", SCchargeHIS, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distHIS')
    # print "HIS", cmd.count_atoms("SC"), SumWSC
    # Sidechain LYS
    nameselect = "/" + molecule + " and resn LYS and name NZ and not resi " + residue
    cmd.select("SC", nameselect)
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distLYS', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeLYS, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC LYS ", str(SClist[i]), " ", SCchargeLYS, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distLYS')
    # print "LYS", cmd.count_atoms("SC"), SumWSC
    # Sidechain MET1
    nameselect = "/" + molecule + " and resn MET and res 1 and not resi " + residue
    cmd.select("SC", nameselect)
    cmd.select("SC", "SC and name N")
    SClist = cmd.identify("SC")
    for i in range(len(SClist)):
        ResDist = cmd.dist(residue + 'distMET1', SGnameselect, molecule + " and id " + str(SClist[i]))
        WSC = fWSC(SCchargeMET1, DieElecSC, ResDist)
        SumWSC = SumWSC + WSC
        if printSC == 'yes':
            print("SC MET1 ", str(SClist[i]), " ", SCchargeMET1, " ", DieElecSC, " ", ResDist, " ", WSC)
    cmd.delete(residue + 'distMET1')
    # print "MET1", cmd.count_atoms("SC"), SumWSC
    cmd.delete("SC")
    return SumWSC


def fWMC(charge, DieElecMC, DistR):
    WMC = 1.39 * math.pow(10, 6) * charge / (DieElecMC * DistR)
    return WMC


def fSumWMCFirst(molecule, SGNameAngle, chain, residue, MCNeighbour, DieElecMC, MCchargeC, MCchargeO, printMC):
    # print "First", MCNeighbour
    SumWMCFirst = 0.0
    SGnameselect = "/" + SGNameAngle + "//" + "/" + "/SG"
    NBnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour)
    cmd.select("MC", NBnameselect)
    MCpdbstr = cmd.get_pdbstr("MC")
    MCsplit = MCpdbstr.split()
    residueName = MCsplit[3]
    # print NBnameselect, residueName
    # Mainchain C
    Cnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/C"
    ResDist = cmd.dist(residue + 'distFirstC', SGnameselect, Cnameselect)
    WMC = fWMC(MCchargeC, DieElecMC, ResDist)
    SumWMCFirst = SumWMCFirst + WMC
    if printMC == 'yes':
        print("MC C ", MCNeighbour, " ", MCchargeC, " ", DieElecMC, " ", ResDist, " ", WMC)
    # Mainchain O
    Onameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/O"
    ResDist = cmd.dist(residue + 'distFirstO', SGnameselect, Onameselect)
    WMC = fWMC(MCchargeO, DieElecMC, ResDist)
    SumWMCFirst = SumWMCFirst + WMC
    if printMC == 'yes':
        print("MC O ", MCNeighbour, " ", MCchargeO, " ", DieElecMC, " ", ResDist, " ", WMC)
    cmd.delete(residue + 'distFirstC')
    cmd.delete(residue + 'distFirstO')
    cmd.delete("MC")
    return SumWMCFirst


def fSumWMCresidue(molecule, SGNameAngle, chain, residue, MCNeighbour, DieElecMC, MCchargeC, MCchargeO, MCchargeN, MCchargeH, AmideName, printMC):
    # print "residue", MCNeighbour
    SumWMCresidue = 0.0
    SGnameselect = "/" + SGNameAngle + "//" + "/" + "/SG"
    NBnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour)
    cmd.select("MC", NBnameselect)
    MCpdbstr = cmd.get_pdbstr("MC")
    MCsplit = MCpdbstr.split()
    residueName = MCsplit[3]
    # print NBnameselect, residueName
    AmideProt = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/H01"
    Hnameselect = "/" + AmideName + "//" + chain + "/" + str(MCNeighbour) + "/H01"
    if cmd.count_atoms(AmideProt) == 0 and cmd.count_atoms(Hnameselect) == 0:
        HbuildSelect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
        cmd.h_add(HbuildSelect)
        cmd.create(AmideName, AmideName + " + " + AmideProt)
        cmd.remove(AmideProt)
    # Mainchain AmideH
    ResDist = cmd.dist(residue + 'distResH', SGnameselect, Hnameselect)
    WMC = fWMC(MCchargeH, DieElecMC, ResDist)
    SumWMCresidue = SumWMCresidue + WMC
    if printMC == 'yes':
        print("MC H ", MCNeighbour, " ", MCchargeH, " ", DieElecMC, " ", ResDist, " ", WMC)
    # Mainchain C
    Cnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/C"
    ResDist = cmd.dist(residue + 'distResC', SGnameselect, Cnameselect)
    WMC = fWMC(MCchargeC, DieElecMC, ResDist)
    SumWMCresidue = SumWMCresidue + WMC
    if printMC == 'yes':
        print("MC C ", MCNeighbour, " ", MCchargeC, " ", DieElecMC, " ", ResDist, " ", WMC)
    # Mainchain O
    Onameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/O"
    ResDist = cmd.dist(residue + 'distResO', SGnameselect, Onameselect)
    WMC = fWMC(MCchargeO, DieElecMC, ResDist)
    SumWMCresidue = SumWMCresidue + WMC
    if printMC == 'yes':
        print("MC O ", MCNeighbour, " ", MCchargeO, " ", DieElecMC, " ", ResDist, " ", WMC)
    # Mainchain N
    Nnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
    ResDist = cmd.dist(residue + 'distResN', SGnameselect, Nnameselect)
    WMC = fWMC(MCchargeN, DieElecMC, ResDist)
    SumWMCresidue = SumWMCresidue + WMC
    if printMC == 'yes':
        print("MC N ", MCNeighbour, " ", MCchargeN, " ", DieElecMC, " ", ResDist, " ", WMC)
    cmd.delete(residue + 'distResH')
    cmd.delete(residue + 'distResC')
    cmd.delete(residue + 'distResO')
    cmd.delete(residue + 'distResN')
    cmd.show("nb_spheres", AmideName)
    cmd.delete("MC")
    return SumWMCresidue


def fSumWMCLast(molecule, SGNameAngle, chain, residue, MCNeighbour, DieElecMC, MCchargeN, MCchargeH, MCchargeProCA, MCchargeProCD, MCchargeProN, AmideName, printMC):
    # print "Last", MCNeighbour
    SumWMCLast = 0.0
    SGnameselect = "/" + SGNameAngle + "//" + "/" + "/SG"
    NBnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour)
    cmd.select("MC", NBnameselect)
    MCpdbstr = cmd.get_pdbstr("MC")
    MCsplit = MCpdbstr.split()
    residueName = MCsplit[3]
    # print NBnameselect, residueName
    if residueName == "PRO":
        # Proline CA
        CAnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/CA"
        ResDist = cmd.dist(residue + 'distLastProCA', SGnameselect, CAnameselect)
        WMC = fWMC(MCchargeProCA, DieElecMC, ResDist)
        SumWMCLast = SumWMCLast + WMC
        if printMC == 'yes':
            print("MC ProCA ", MCNeighbour, " ", MCchargeProCA, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Proline CD
        CDnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/CD"
        ResDist = cmd.dist(residue + 'distLastProCD', SGnameselect, CDnameselect)
        WMC = fWMC(MCchargeProCD, DieElecMC, ResDist)
        SumWMCLast = SumWMCLast + WMC
        if printMC == 'yes':
            print("MC ProCD ", MCNeighbour, " ", MCchargeProCD, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Proline N
        Nnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
        ResDist = cmd.dist(residue + 'distLastProN', SGnameselect, Nnameselect)
        WMC = fWMC(MCchargeProN, DieElecMC, ResDist)
        SumWMCLast = SumWMCLast + WMC
        if printMC == 'yes':
            print("MC ProN ", MCNeighbour, " ", MCchargeProN, " ", DieElecMC, " ", ResDist, " ", WMC)
    else:
        AmideProt = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/H01"
        Hnameselect = "/" + AmideName + "//" + chain + "/" + str(MCNeighbour) + "/H01"
        if cmd.count_atoms(AmideProt) == 0 and cmd.count_atoms(Hnameselect) == 0:
            HbuildSelect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
            cmd.h_add(HbuildSelect)
            cmd.create(AmideName, AmideName + " + " + AmideProt)
            cmd.remove(AmideProt)
        # Mainchain AmideH
        ResDist = cmd.dist(residue + 'distLastH', SGnameselect, Hnameselect)
        WMC = fWMC(MCchargeH, DieElecMC, ResDist)
        SumWMCLast = SumWMCLast + WMC
        if printMC == 'yes':
            print("MC H ", MCNeighbour, " ", MCchargeH, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Mainchain N
        Nnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
        ResDist = cmd.dist(residue + 'distLastN', SGnameselect, Nnameselect)
        WMC = fWMC(MCchargeN, DieElecMC, ResDist)
        SumWMCLast = SumWMCLast + WMC
        if printMC == 'yes':
            print("MC N ", MCNeighbour, " ", MCchargeN, " ", DieElecMC, " ", ResDist, " ", WMC)
    cmd.delete(residue + 'distLastProCA')
    cmd.delete(residue + 'distLastProCD')
    cmd.delete(residue + 'distLastProN')
    cmd.delete(residue + 'distLastH')
    cmd.delete(residue + 'distLastN')
    cmd.show("nb_spheres", AmideName)
    cmd.delete("MC")
    return SumWMCLast


def fSumWMC(molecule, SGNameAngle, chain, residue, MCNeighbour, DieElecMC, MCchargeC, MCchargeO, MCchargeN, MCchargeH, MCchargeProCA, MCchargeProCD, MCchargeProN, AmideName, printMC):
    # print "chain", MCNeighbour
    SumWMC = 0.0
    SGnameselect = "/" + SGNameAngle + "//" + "/" + "/SG"
    NBnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour)
    cmd.select("MC", NBnameselect)
    MCpdbstr = cmd.get_pdbstr("MC")
    MCsplit = MCpdbstr.split()
    residueName = MCsplit[3]
    # print NBnameselect, residueName
    if residueName == "PRO":
        # Proline CA
        CAnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/CA"
        ResDist = cmd.dist(residue + 'distProCA', SGnameselect, CAnameselect)
        WMC = fWMC(MCchargeProCA, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC ProCA ", MCNeighbour, " ", MCchargeProCA, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Proline CD
        CDnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/CD"
        ResDist = cmd.dist(residue + 'distProCD', SGnameselect, CDnameselect)
        WMC = fWMC(MCchargeProCD, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC ProCD ", MCNeighbour, " ", MCchargeProCD, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Proline N
        Nnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
        ResDist = cmd.dist(residue + 'distProN', SGnameselect, Nnameselect)
        WMC = fWMC(MCchargeProN, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC ProN ", MCNeighbour, " ", MCchargeProN, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Proline C
        Cnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/C"
        ResDist = cmd.dist(residue + 'distProC', SGnameselect, Cnameselect)
        WMC = fWMC(MCchargeC, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC ProC ", MCNeighbour, " ", MCchargeC, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Proline O
        Onameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/O"
        ResDist = cmd.dist(residue + 'distProO', SGnameselect, Onameselect)
        WMC = fWMC(MCchargeO, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC ProO ", MCNeighbour, " ", MCchargeO, " ", DieElecMC, " ", ResDist, " ", WMC)
    else:
        AmideProt = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/H01"
        Hnameselect = "/" + AmideName + "//" + chain + "/" + str(MCNeighbour) + "/H01"
        if cmd.count_atoms(AmideProt) == 0 and cmd.count_atoms(Hnameselect) == 0:
            HbuildSelect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
            cmd.h_add(HbuildSelect)
            cmd.create(AmideName, AmideName + " + " + AmideProt)
            cmd.remove(AmideProt)
        # Mainchain AmideH
        ResDist = cmd.dist(residue + 'distH', SGnameselect, Hnameselect)
        WMC = fWMC(MCchargeH, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC H ", MCNeighbour, " ", MCchargeH, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Mainchain C
        Cnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/C"
        ResDist = cmd.dist(residue + 'distC', SGnameselect, Cnameselect)
        WMC = fWMC(MCchargeC, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC C ", MCNeighbour, " ", MCchargeC, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Mainchain O
        Onameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/O"
        ResDist = cmd.dist(residue + 'distO', SGnameselect, Onameselect)
        WMC = fWMC(MCchargeO, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC O ", MCNeighbour, " ", MCchargeO, " ", DieElecMC, " ", ResDist, " ", WMC)
        # Mainchain N
        Nnameselect = "/" + molecule + "//" + chain + "/" + str(MCNeighbour) + "/N"
        ResDist = cmd.dist(residue + 'distN', SGnameselect, Nnameselect)
        WMC = fWMC(MCchargeN, DieElecMC, ResDist)
        SumWMC = SumWMC + WMC
        if printMC == 'yes':
            print("MC N ", MCNeighbour, " ", MCchargeN, " ", DieElecMC, " ", ResDist, " ", WMC)
    cmd.delete(residue + 'distProCA')
    cmd.delete(residue + 'distProCD')
    cmd.delete(residue + 'distProN')
    cmd.delete(residue + 'distProC')
    cmd.delete(residue + 'distProO')
    cmd.delete(residue + 'distH')
    cmd.delete(residue + 'distC')
    cmd.delete(residue + 'distO')
    cmd.delete(residue + 'distN')
    cmd.show("nb_spheres", AmideName)
    cmd.delete("MC")
    return SumWMC


def fBoltzSingleState(SumMCSC, AvogadroR, Temp):
    BoltzSingleState = math.exp(-SumMCSC / (AvogadroR * Temp))
    return BoltzSingleState


def fpKm1(DeltapKMCSC, pK1):
    pKm1 = DeltapKMCSC + pK1
    return pKm1


def fpKm2(DeltapKMCSC, DeltapKB, pK2):
    pKm2 = DeltapKMCSC + DeltapKB + pK2
    return pKm2


def fFracCys(pKm, pH):
    FracCys = 1.0 / (math.pow(10, (pKm - pH)) + 1)
    return FracCys


def AtomVector(AtomStart, AtomEnd):
    PosStart = cmd.get_atom_coords(AtomStart)
    PosEnd = cmd.get_atom_coords(AtomEnd)
    VectorDiff = [(PosEnd[0] - PosStart[0]), (PosEnd[1] - PosStart[1]), (PosEnd[2] - PosStart[2])]
    return VectorDiff


def pairfitCys(Cysmolecule, molecule, chain, residue):
    RN = "/" + Cysmolecule + "//" + "/" + "/N"
    PN = "/" + molecule + "//" + chain + "/" + residue + "/N"
    RCA = "/" + Cysmolecule + "//" + "/" + "/CA"
    PCA = "/" + molecule + "//" + chain + "/" + residue + "/CA"
    RC = "/" + Cysmolecule + "//" + "/" + "/C"
    PC = "/" + molecule + "//" + chain + "/" + residue + "/C"
    RCB = "/" + Cysmolecule + "//" + "/" + "/CB"
    PCB = "/" + molecule + "//" + chain + "/" + residue + "/CB"
    cmd.select("CBatom", PCB)
    CBatomNr = cmd.count_atoms("CBatom")
    # If PRO or GLY, then only fit N, CA, C atoms
    if CBatomNr == 0:
        cmd.pair_fit(RN, PN, RCA, PCA, RC, PC)
    else:
        # cmd.pair_fit(RN,PN,RCA,PCA,RC,PC,RCB,PCB)
        cmd.pair_fit(RN, PN, RCA, PCA, RC, PC)
    cmd.delete("CBatom")


def CheckDimer(dihedSG, molecule, chain, residue):
    breakDimer = "no"
    nameselect = "(/" + molecule + "//" + chain + " and name SG and not /" + molecule + "//" + chain + "/" + residue + ") within 5 of " + dihedSG
    cmd.select(str(molecule) + str(residue) + "Dimer", nameselect)
    DimerSG = cmd.count_atoms(str(molecule) + str(residue) + "Dimer")
    if DimerSG > 0:
        print("####################################################")
        print("########### WARNING: SG in near detected ###########")
        print("########### Is this a dimer?             ###########")
        print("####################################################")
        cmd.select(str(molecule) + str(residue) + "Dimer", "byres " + str(molecule) + str(residue) + "Dimer")
        cmd.show("sticks", str(molecule) + str(residue) + "Dimer")
        breakDimer = "yes"
    else:
        cmd.delete(str(molecule) + str(residue) + "Dimer")
    return breakDimer
