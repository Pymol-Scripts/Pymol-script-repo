'''
Described at PyMOL wiki:
http://www.pymolwiki.org/index.php/rotkit

-------------------------------------------------------------------------------
 Name:		rotkit.py   examples
 Purpose:      To rotate molecules easier

 Author:      Troels Linnet

 Created:     30/08/2011
 Copyright:   (c) Troels Linnet 2011
 Licence:     Free
-------------------------------------------------------------------------------
'''

from __future__ import print_function

from pymol import cmd
import math
import os
import platform


def printMat(matrix):
    print(("%s %s %s %s \n%s %s %s %s \n%s %s %s %s \n%s %s %s %s" % (matrix[0], matrix[1], matrix[2], matrix[3], matrix[4], matrix[5], matrix[6], matrix[7], matrix[8], matrix[9], matrix[10], matrix[11], matrix[12], matrix[13], matrix[14], matrix[15])))
    return None


def getxyz(Sel):
    if type(Sel) == list and len(Sel) == 3:
        return Sel, "listXYZ"
    if type(Sel) == str and Sel[0] == "[" and Sel[-1] == "]":
        Selsplit = list(Sel[1:-1].split(","))
        Selsplit = [float(x) for x in Selsplit]
        return Selsplit, "strXYZ"
    if type(Sel) == str:
        pos = cmd.get_atom_coords(Sel)
        return pos, "selXYZ"


def vector(Sel1, Sel2):
    PosSel1 = getxyz(Sel1)[0]
    PosSel2 = getxyz(Sel2)[0]
    vectorcalc = [PosSel2[0] - PosSel1[0], PosSel2[1] - PosSel1[1], PosSel2[2] - PosSel1[2]]
    return(vectorcalc)


def vectoradd(Sel1, Sel2):
    PosSel1 = getxyz(Sel1)[0]
    PosSel2 = getxyz(Sel2)[0]
    vectorcalc = [PosSel1[0] + PosSel2[0], PosSel1[1] + PosSel2[1], PosSel1[2] + PosSel2[2]]
    return(vectorcalc)


def vectorstr(vector):
    return("[%s,%s,%s]" % (vector[0], vector[1], vector[2]))


def transmat(vector, dist=1):
    mat = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, dist * vector[0], dist * vector[1], dist * vector[2], 1]
    return(mat)


def unitvector(vector):
    vectorlen = math.sqrt(math.pow(vector[0], 2) + math.pow(vector[1], 2) + math.pow(vector[2], 2))
    vectordiv = [vector[0] / vectorlen, vector[1] / vectorlen, vector[2] / vectorlen]
    return(vectordiv, vectorlen)


def radangle(angle):
    return(math.radians(angle))


def rotmat(angle, vectornorm, pointcoord):
    # From: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ Section 6.2
    u, v, w = vectornorm
    a, b, c = pointcoord
    makerotmat = [(math.pow(u, 2) + (math.pow(v, 2) + math.pow(w, 2)) * math.cos(angle)),
                  (u * v * (1 - math.cos(angle)) - w * math.sin(angle)),
                  (u * w * (1 - math.cos(angle)) + v * math.sin(angle)),
                  ((a * (math.pow(v, 2) + math.pow(w, 2)) - u * (b * v + c * w)) * (1 - math.cos(angle)) + (b * w - c * v) * math.sin(angle)),
                  (u * v * (1 - math.cos(angle)) + w * math.sin(angle)),
                  (math.pow(v, 2) + (math.pow(u, 2) + math.pow(w, 2)) * math.cos(angle)),
                  (v * w * (1 - math.cos(angle)) - u * math.sin(angle)),
                  ((b * (math.pow(u, 2) + math.pow(w, 2)) - v * (a * u + c * w)) * (1 - math.cos(angle)) + (c * u - a * w) * math.sin(angle)),
                  (u * w * (1 - math.cos(angle)) - v * math.sin(angle)),
                  (v * w * (1 - math.cos(angle)) + u * math.sin(angle)),
                  (math.pow(w, 2) + (math.pow(u, 2) + math.pow(v, 2)) * math.cos(angle)),
                  ((c * (math.pow(u, 2) + math.pow(v, 2)) - w * (a * u + b * v)) * (1 - math.cos(angle)) + (a * v - b * u) * math.sin(angle)),
                  (0), (0), (0), (1), ]
    return(makerotmat)


def rotateline(Pos1, Pos2, degangle, molecule):
    diffvector = vector(Pos1, Pos2)
    uvector = unitvector(diffvector)[0]
    xyz = getxyz(Pos2)[0]
    rmat = rotmat(radangle(float(degangle)), uvector, xyz)
    cmd.transform_selection(molecule, rmat)
    return(None)
cmd.extend("rotateline", rotateline)


def mutate(molecule, chain, resi, target="CYS", mutframe="1"):
    target = target.upper()
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    cmd.get_wizard().set_mode("%s" % target)
    selection = "/%s//%s/%s" % (molecule, chain, resi)
    cmd.get_wizard().do_select(selection)
    cmd.frame(str(mutframe))
    cmd.get_wizard().apply()
    # cmd.set_wizard("done")
    cmd.set_wizard()
    # cmd.refresh()
#### Example in pymol
# python
#MutList = [["5NT","A",308],["5NT","A",513],["5NT","B",513]]
# for p,c,r in MutList:
#    rotkit.mutate(p, chain=c, resi=r, target="CYS", mutframe=1)
# Have do mutate first before selecting, or else it only select the lase?
# for p,c,r in MutList:
#    cmd.select("%s%s%s"%(p,c,r),"/%s//%s/%s"%((p,c,r)))
# python end
cmd.extend("mutate", mutate)


def toline(Pos1, Pos2, atom, molecule, dist=1):
    dist = float(dist)
    diffvector = vector(atom, Pos2)
    move = transmat(diffvector)
    cmd.transform_selection("%s" % molecule, move)
    diffvector = vector(Pos1, Pos2)
    uvector = unitvector(diffvector)[0]
    move = transmat(uvector, dist)
    cmd.transform_selection("%s" % molecule, move)
    return(None)
cmd.extend("toline", toline)


def crossprod(Vector1, Vector2):
    return([Vector1[1] * Vector2[2] - Vector1[2] * Vector2[1], Vector1[2] * Vector2[0] - Vector1[0] * Vector2[2], Vector1[0] * Vector2[1] - Vector1[1] * Vector2[0]])


def crosspoint(Pos1, crossprod):
    Imp1 = getxyz(Pos1)[0]
    Imp2 = getxyz(crossprod)[0]
    return([Imp1[0] + Imp2[0], Imp1[1] + Imp2[1], Imp1[2] + Imp2[2]])


def VectorToMatrix(Vector, MatColRank=4):
    try:
        import numpy
    except ImportError:
        from modules import numpy
    nextrow = list(range(MatColRank, MatColRank ** 2, MatColRank))
    rowsall = []
    rowcurrent = []
    for i in range(len(Vector)):
        if i in nextrow:
            rowsall.append(rowcurrent)
            rowcurrent = []
            rowcurrent.append(Vector[i])
        else:
            rowcurrent.append(Vector[i])
    print(rowsall)
    return(numpy.matrix(rowsall))


def findMinMax(datalist, index):
    minimum = datalist[0][index]
    maximum = datalist[0][index]
    datacolumn = []
    for l in datalist:
        datacolumn.append(l[index])
        if l[index] < minimum:
            minimum = l[index]
        if l[index] > maximum:
            maximum = l[index]
    return(minimum, maximum, datacolumn)


def createdirs(dirname):
    if platform.system() == 'Windows':
        Newdir = os.getcwd() + "\\%s\\" % dirname
    if platform.system() == 'Linux':
        Newdir = os.getcwd() + "/%s/" % dirname
    if not os.path.exists(Newdir):
        os.makedirs(Newdir)
    return(Newdir)


def makehistogram(datalist, dataname="Histogram", datalistindex=2, nrbins=100, binrange=[0, 0]):
    try:
        import numpy
    except ImportError:
        from modules import numpy
    fileout_name = "%s" % dataname + ".dat"
    fileout_write = open(fileout_name, "w")
    gnuplot_write = open("%s" % dataname + ".plt", "w")
    datacolumnMin, datacolumnMax, datacolumn = findMinMax(datalist, datalistindex)
    if binrange[1] == 0:
        xdelta = datacolumnMax - datacolumnMin
        binrange[0] = datacolumnMin - xdelta / 2
        binrange[1] = datacolumnMax + xdelta / 2
    # print binrange
    (n, binval) = numpy.histogram(datacolumn, bins=int(nrbins), range=(binrange[0], binrange[1]), normed=False)
    binwidthDist = (binrange[1] - binrange[0]) / nrbins
    DistHist = []
    # Normalize the histogram
    for i in range(len(n)):
        DistHist.append([binval[i], n[i], float(n[i]) / len(datacolumn)])
    # print DistHist
    DistHistMin, DistHistMax, tmp = findMinMax(DistHist, 2)
    # print DistHistMin,DistHistMax
    ydelta = DistHistMax - DistHistMin
    # Now write the output
    fileout_write.write("#Datapoints=%s" % len(datacolumn) + "\n")
    fileout_write.write("#Dist[Ang]    Frequency[#]    Probability" + "\n")
    for dp in DistHist:
        textline = "%4.2f    %5i    %18.5f" % (dp[0], dp[1], dp[2])
        fileout_write.write(textline + "\n")
    gnuplot_write.write('cd "%s"' % os.getcwd() + '\n')
    gnuplot_write.write('set term postscript eps enhanced color' + '\n')
    gnuplot_write.write('' + '\n')
    gnuplot_write.write('set style line 1 lt 1 lw 3 linecolor rgb "red"' + '\n')
    gnuplot_write.write('set style line 2 lt 1 lw 3 linecolor rgb "green"' + '\n')
    gnuplot_write.write('set style line 3 lt 1 lw 3 linecolor rgb "blue"' + '\n')
    gnuplot_write.write('set style line 4 lt 1 lw 3 linecolor rgb "pink"' + '\n')
    gnuplot_write.write('set style line 5 lt 1 lw 0 linecolor rgb "red"' + '\n')
    gnuplot_write.write('binwidthDist=%s' % binwidthDist + '\n')
    gnuplot_write.write('set title "Normalized distance histogram"' + '\n')
    gnuplot_write.write('set xlabel "Distance [Ang]"' + '\n')
    gnuplot_write.write('set xrange[%s:%s]' % (binrange[0], binrange[1]) + '\n')
    gnuplot_write.write('set ylabel "Density"' + '\n')
    gnuplot_write.write('set yrange[%s:%s]' % (DistHistMin, DistHistMax / binwidthDist) + '\n')
    gnuplot_write.write('set ytics nomirror' + '\n')
    gnuplot_write.write('set y2label "Integrated Bin probability"' + '\n')
    gnuplot_write.write('set y2range[%s:%s]' % (DistHistMin, DistHistMax) + '\n')
    gnuplot_write.write('set y2tics border' + '\n')
    gnuplot_write.write('' + '\n')
    gnuplot_write.write('A=1' + '\n')
    gnuplot_write.write('sigma2=%s' % (xdelta / 6.0) + '\n')
    gnuplot_write.write('center=%s' % (datacolumnMin + xdelta / 2) + '\n')
    gnuplot_write.write('g(x) = (A*1.0/sqrt(2*pi*sigma2))*exp(-(x-center)**2/(2*sigma2))' + '\n')
    gnuplot_write.write('fit g(x) "%s" using 1:($3/binwidthDist) via sigma2,center' % fileout_name + '\n')
    gnuplot_write.write('' + '\n')
    gnuplot_write.write('set label "g(x)= A*1/(sqrt(2{/Symbol p}{/Symbol s}^2)) * exp(-(x-{/Symbol m})^2/(2{/Symbol s}^2))" at graph 0.05, graph 0.85' + '\n')
    gnuplot_write.write('set label "A= %g", A at graph 0.05, graph 0.80' + '\n')
    gnuplot_write.write('set label "{/Symbol s}= %g", sqrt(sigma2) at graph 0.05, graph 0.75' + '\n')
    gnuplot_write.write('set label "{/Symbol m}= %g", center at graph 0.05, graph 0.70' + '\n')
    gnuplot_write.write('' + '\n')
    gnuplot_write.write('set label "Binwidth= %g", binwidthDist at graph 0.05, graph 0.95' + '\n')
    gnuplot_write.write('set label "Datapoints= %g", ' + '%s' % len(datacolumn) + ' at graph 0.05, graph 0.90' + '\n')
    gnuplot_write.write('set output "%s.eps"' % (fileout_name[:-4]) + '\n')
    gnuplot_write.write('plot "%s" using 1:($3/binwidthDist) title "%s" with boxes fs solid 0.4 noborder,\\' % (fileout_name, fileout_name[:-4]) + '\n')
    gnuplot_write.write('g(x) title "Fitted normal distribution g(x)" lw 4,\\' + '\n')
    gnuplot_write.write('"%s" using 1:3 title "" with histeps ls 5 axis x1y2' % fileout_name + '\n')

    fileout_write.close()
    gnuplot_write.close()
    return(DistHist)
