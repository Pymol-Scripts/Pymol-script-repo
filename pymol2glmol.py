'''
Pymol to GLmol exporter

Written by biochem_fan, 2011
'''

from pymol import cmd
from math import cos, sin, pi, sqrt, acos, asin, atan2
import os

def compactSeq(seq):
    seq.sort()
    ret = []
    prev = -9999
    start = -9999
    seq.append(-1)
    i = 0
    while (i < len(seq) - 1):
        if (start >= 0):
            if (seq[i] + 1 != seq[i + 1]):
                ret.append("%d-%d" % (start, seq[i]))
                start = -999
        else:
            if (seq[i] + 1 != seq[i + 1]):
                start = -999
                ret.append(str(seq[i]))
            else:
                start = seq[i]
        i += 1
    return ','.join(ret)

def parseObjMol(obj):
    name = obj[0]
    ids = []
    sphere = []
    trace = []
    ribbon = []
    stick = []
    surface = []
    line = []
    cross = []
    smallSphere = []
    helix = []
    sheet = []
    colors = {}
    for atom in obj[5][7]:
        rep = atom[20]
        serial = atom[22]
        ss = atom[10]
        bonded = (atom[25] == 1)
        if (rep[5] == 1):
            ribbon.append(serial)
        if (rep[1] == 1):
            sphere.append(serial)
        if (rep[2] == 1):
            surface.append(serial)
        if (rep[7] == 1):
            line.append(serial)
        if (rep[6] == 1):
            trace.append(serial)
        if (rep[4] == 1 and not bonded):
            smallSphere.append(serial)
        if (rep[11] == 1 and not bonded):
            cross.append(serial)
        if (rep[0] == 1 and bonded):
            stick.append(serial)
        if (ss == 'S'):
            sheet.append(serial)
        if (ss == 'H'):
            helix.append(serial)

        c =  cmd.get_color_tuple(atom[21])
        if (not c in colors):
            colors[c] = []
        colors[c].append(serial)
        ids.append("ID %d is %s in resi %s %s at chain %s"\
                       % (atom[22], atom[6], atom[3], atom[5], atom[1]))

    for c in colors.iterkeys(): # TODO: better compression
        colors[c] = compactSeq(colors[c])

    ret = ''
    ret += "\nsheet:" + compactSeq(sheet)
    ret += "\nhelix:" + compactSeq(helix)
    ret += "\nsurface:" + compactSeq(surface)
    ret += "\nsphere:" + compactSeq(sphere)
    ret += "\ntrace:" + compactSeq(trace)
    ret += "\nribbon:" + compactSeq(ribbon)
    ret += "\nstick:" + compactSeq(stick)
    ret += "\nline:" + compactSeq(line)
    ret += "\nsmallSphere:" + compactSeq(smallSphere)
    ret += "\ncross:" + compactSeq(cross)
    for c in colors.iterkeys():
        ret += "\ncolor:%.3f,%.3f,%.3f:%s" % (c[0], c[1], c[2], colors[c])
    return ret

def parseDistObj(obj):
    if (obj[5][0][3][10] != 1): # 'show dashed' flag
        return ""
    N = obj[5][2][0][0]
    points = obj[5][2][0][1]
    ret = []
    for p in points:
        ret.append("%.3f" % p)
    color = cmd.get_color_tuple(obj[5][0][2]);
    return "\ndists:%.3f,%.3f,%.3f:" % color + ','.join(ret)

def dump_rep(name):
    if 'PYMOL_GIT_MOD' in os.environ:
        import shutil
        try:
            shutil.copytree(os.path.join(os.environ['PYMOL_GIT_MOD'],'pymol2glmol','js'),os.path.join(os.getcwd(),'js'))
        except OSError:
            pass

    names = cmd.get_session()['names']
    cmd.set('pdb_retain_ids', 1)

    ret = ''
    for obj in names:
        if (obj == None):
            continue
        if (obj[2] == 0): # not visible
            continue
        if (obj[1] == 0 and obj[4] == 1 and obj[0] == name):
            ret += parseObjMol(obj)
        if (obj[1] == 0 and obj[4] == 4): # currently all dist objects are exported
            ret += parseDistObj(obj)

    cmd.turn('z', 180)
    view = cmd.get_view()
    cmd.turn('z', 180)
    cx = -view[12]; cy = -view[13]; cz = -view[14]
    cameraZ = - view[11] - 150;
    fov = float(cmd.get("field_of_view"))
    fogStart = float(cmd.get("fog_start"))
    slabNear = view[15] + view[11]
    slabFar = view[16] + view[11]
    ret += "\nview:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f" % \
        (cx, cy, cz, cameraZ, slabNear, slabFar, fogStart, fov)
    for i in range(9):
        ret += ",%.3f" % view[i]

    bgcolor = cmd.get_setting_tuple('bg_rgb')[1]
    ret += "\nbgcolor:%02x%02x%02x" % (int(255 * float(bgcolor[0])), \
              int(255 * float(bgcolor[1])), int(255 * float(bgcolor[2])))
    if 'PYMOL_GIT_MOD' in os.environ:
        template = open(os.path.join(os.environ['PYMOL_GIT_MOD'],'pymol2glmol','imported.html')).read().\
            replace("###INCLUDE_PDB_FILE_HERE###", cmd.get_pdbstr(name)).\
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)
    else:
        template = open('imported.html').read().\
            replace("###INCLUDE_PDB_FILE_HERE###", cmd.get_pdbstr(name)).\
            replace('###INCLUDE_REPRESENTATION_HERE###', ret)
        
    f = open(name + '.html', 'w')
    f.write(template)
    f.close()

cmd.extend('pymol2glmol', dump_rep)
cmd.auto_arg[0]['pymol2glmol'] = [cmd.object_sc, 'object', '']
