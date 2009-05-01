##################################################################################
#movs.py - Math and animation routines for slerpy
#Copyright (C) 2006 Joel Bard
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#################################################################################
 
from pymol import cmd,stored
from math import *
from string import split
 
def rmat2quat( M ):
    #M is a list of 9 values being the elements of the rotation matrix in row order
    T = M[0] + M[4] + M[8] + 1
    print "Trace ",T
    if T>0:
        S = 0.5 / sqrt(T)
        W = 0.25/S
        X = (M[7] - M[5])*S
        Y = (M[2] - M[6])*S
        Z = (M[3] - M[1])*S
    else:
        if( M[0] > M[4] and M[0] > M[8] ):
            S = sqrt( 1.0 + M[0] - M[4] - M[8]) * 2
            X = 0.25 * S
            Y = (M[1] + M[3])/S
            Z = (M[2] + M[6])/S
            W = (M[5] - M[7])/S
        elif M[4] > M[8]:
            S = sqrt( 1.0 + M[4] - M[0] - M[8] ) * 2
            X = (M[1] + M[3])/S
            Y = 0.25 * S
            Z = (M[5] + M[7])/S
            W = (M[2] - M[6])/S
        else:
            S = sqrt( 1.0 + M[8] - M[0] - M[4]) * 2
            X = (M[2] + M[6])/S
            Y = (M[5] + M[7])/S
            Z = 0.25 * S
            W = (M[1] - M[3])/S
    return [X,Y,Z,W]
 
def quat2rmat( Q ):
    #Q is a list of 4 values being the quaternion X Y Z W
    X=Q[0]
    Y=Q[1]
    Z=Q[2]
    W=Q[3]
    xx = X*X
    xy = X*Y
    xz = X*Z
    xw = X*W
    yy = Y*Y
    yz = Y*Z
    yw = Y*W
    zz = Z*Z
    zw = Z*W
 
    M= [1.0]*9
    M[0] = 1 - 2 * ( yy + zz )
    M[1] = 2 * ( xy - zw )
    M[2] = 2 * ( xz + yw )
    M[3] = 2 * ( xy + zw )
    M[4] = 1 - 2 * ( xx + zz )
    M[5] = 2 * ( yz - xw )
    M[6] = 2 * ( xz - yw )
    M[7] = 2 * ( yz + xw )
    M[8] = 1 - 2 * ( xx + yy )
    return M
 
def quatconj( Q ):
    return [-Q[0],-Q[1],-Q[2],Q[3]]
 
def quatmag( Q ):
    s = 0.0
    QC = quatconj(Q)
    for x in range(0,4):
        s += Q[x]*Q[x]
    print s
    return sqrt(s)
 
def quatnorm( Q ):
    m = quatmag( Q )
    return [q/m for q in Q]
 
def quatdotprod( q1, q2 ):
    dp = 0
    for i in range(0,4):
        dp += q1[i]*q2[i]
    return dp
 
def vectnorm( V ):
    mag = 0.0
    for x in V:
        mag += x*x
    mag = sqrt(mag)
    return [x/mag for x in V]
 
def quat2axisangle( Q ):
    #returns list where 0..2 are rot axis and 3 is angle
    qn = quatnorm( Q )
    cos_a = Q[3]
    angle = acos( cos_a ) * 2
    sin_a = sqrt( 1.0 - cos_a * cos_a )
    if( fabs( sin_a ) < 0.000005 ):
        sin_a = 1
    ax_an = [ q/sin_a for q in Q[0:3] ]
    ax_an.append( angle )
    return ax_an
 
def axisangle2quat( ax_an ):
    #ax_an is a list with axis coordinates followed by rotation angle
    axn = vectnorm( ax_an[:3] )
    angle = ax_an[3]
    sin_a = sin( angle / 2 )
    cos_a = cos( angle / 2 )
    Q = [ x * sin_a for x in axn ]
    Q.append( cos_a )
    return Q
 
def rmat2axisangle( M ):
    q = rmat2quat( M )
    return quat2axisangle( q )
 
def axisangle2rmat( a ):
    q = axisangle2quat( a )
    return quat2rmat( q )
 
def animate_transition( start_view, end_view, nframes, first_frame, settings = [] ):
    #print "Views"
    #print start_view,'\n',end_view
 
    cview = start_view[:]
    cmd.set_view( start_view )
 
    #get initial and final quaternions for interpolation
    #print start_view[0:9]
    #get quaternions
    qstart = rmat2quat( start_view[0:9] )
    qend = rmat2quat( end_view[0:9] )
 
    #test for long way vs. short way
    if( quatdotprod( qstart,qend ) < 0 ):
        qstart = [-q for q in qstart]
 
    axan_start = quat2axisangle( qstart )
    axan_end = quat2axisangle( qend )
 
    axan_cur = axan_start[:]
    frame_start = first_frame
    frame_end = frame_start + nframes
    doFade = 0
    doSetting = 0
    if len( settings ) == 4:
        settingName, selection, startVal, endVal = settings
        settingStep = (endVal-startVal)/float(nframes)
        print "Setting step ", settingStep
        doSetting = 1
    elif len( settings ) == 3:
        startVisSelection, endVisSelection, sticksOnly = settings
        settingStep = 1.0/float(nframes)
        doFade = 1
    for f in range( frame_start , frame_end):
        #get rotmat
        #using angle axis
 
        for i in range(0,4):
            axan_cur[i] = axan_cur[i] + (axan_end[i]-axan_start[i])/nframes
        newmat = axisangle2rmat( axan_cur )
        #print cview
        for i in range(0,9):
            cview[i] = newmat[i]
 
        mdo_cmd = "set_view (["
        for i in range(0,18):
            if( i>8 ):
                cview[i] = cview[i]+(end_view[i]-start_view[i])/nframes
            mdo_cmd += "%12.7f,"% cview[i]
        mdo_cmd = mdo_cmd[:-1]+"])"
        if doSetting:       
            val = float(f-frame_start)*settingStep + startVal 
            print val;
            mdo_cmd += "; set %s, %f, %s" % (settingName, val, selection)
            print mdo_cmd;
        #print "mdo ", mdo_cmd
        if doFade:
            val = float(f-frame_start)*settingStep
            otherVal = 1.0 - val
            mdo_cmd += "; set stick_transparency, %f, %s; set stick_transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
            if not sticksOnly:
                #comment out surface transparency interpolation due to problem with transparent sticks in front of 
                #transparent surfaces (get black holes)
               # mdo_cmd += "; set transparency, %f, %s; set transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
                mdo_cmd += "; set cartoon_transparency, %f, %s; set cartoon_transparency, %f, %s" % ( val, startVisSelection, otherVal, endVisSelection )
        cmd.mdo(f,mdo_cmd)
