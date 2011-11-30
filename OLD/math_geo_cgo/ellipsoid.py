from pymol.cgo import BEGIN, COLOR, TRIANGLES, VERTEX, NORMAL, END
from pymol import cmd
 
def signOfFloat(f):
        if f < 0: return -1
        if f > 0: return 1
        return 0
 
def sqC(v, n):
        return signOfFloat(math.cos(v)) *  math.pow(math.fabs(math.cos(v)), n)
 
def sqCT(v, n, alpha):
        return alpha + sqC(v, n)
 
def sqS(v, n):
        return signOfFloat(math.sin(v)) * math.pow(math.fabs(math.sin(v)), n)
 
def sqEllipsoid(x, y, z, a1, a2, a3, u, v, n, e):
        x = a1 * sqC(u, n) * sqC(v, e) + x
        y = a2 * sqC(u, n) * sqS(v, e) + y
        z = a3 * sqS(u, n) + z
        nx = sqC(u, 2 - n) * sqC(v, 2 - e) / a1
        ny = sqC(u, 2 - n) * sqS(v, 2 - e) / a2
        nz = sqS(u, 2 - n) / a3
        return x, y, z, nx, ny, nz
 
def sqToroid(x, y, z, a1, a2, a3, u, v, n, e, alpha):
        a1prime = 1.0 / (a1 + alpha)
        a2prime = 1.0 / (a2 + alpha)
        a3prime = 1.0 / (a3 + alpha)
        x = a1prime * sqCT(u, e, alpha) * sqC(v, n)
        y = a2prime * sqCT(u, e, alpha) * sqS(v, n)
        z = a3prime * sqS(u, e)
        nx = sqC(u, 2 - e) * sqC(v, 2 - n) / a1prime
        ny = sqC(u, 2 - e) * sqS(v, 2 - n) / a2prime
        nz = sqS(u, 2 - e) / a3prime
        return x, y, z, nx, ny, nz
 
def makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, n, e, u1, u2, v1, v2, u_segs, v_segs, color=[0.5, 0.5, 0.5]):
 
        r, g, b = color
 
        # Calculate delta variables */
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs
 
        o = [ BEGIN, TRIANGLES ]
 
        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop */
                V = v1
                for X in range(0, v_segs):
                        # VERTEX #1 */
                        x1, y1, z1, n1x, n1y, n1z = sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, e)
                        x2, y2, z2, n2x, n2y, n2z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, e)
                        x3, y3, z3, n3x, n3y, n3z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e)
                        x4, y4, z4, n4x, n4y, n4z = sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, e)
 
                        o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
 
                        # Update variables for next loop */
                        V += dV
                # Update variables for next loop */
                U += dU
        o.append(END)
        return o
 
def makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, n, e, u1, u2, v1, v2, u_segs, v_segs, color=[0.5, 0.5, 0.5]):
 
        r, g, b = color
 
        # Calculate delta variables */
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs
 
        o = [ BEGIN, TRIANGLES ]
 
        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop */
                V = v1
                for X in range(0, v_segs):
                        # VERTEX #1 */
                        x1, y1, z1, n1x, n1y, n1z = sqToroid(x, y, z, a1, a2, a3, U, V, n, e, alpha)
                        x2, y2, z2, n2x, n2y, n2z = sqToroid(x, y, z, a1, a2, a3, U + dU, V, n, e, alpha)
                        x3, y3, z3, n3x, n3y, n3z = sqToroid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e, alpha)
                        x4, y4, z4, n4x, n4y, n4z = sqToroid(x, y, z, a1, a2, a3, U, V + dV, n, e, alpha)
 
                        o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
 
                        # Update variables for next loop */
                        V += dV
                # Update variables for next loop */
                U += dU
        o.append(END)
        return o
 
def makeEllipsoid(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
 
def makeCylinder(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
 
def makeSpindle(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 1.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
 
def makeDoublePyramid(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 2.0, 2.0, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
 
def makePillow(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 1.0, 0.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)
 
def makeRoundCube(x, y, z, a1, a2, a3):
                return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, 0.2, 0.2, -math.pi / 2, math.pi / 2, -math.pi, math.pi, 10, 10)
 
def makeToroid(x, y, z, a1, a2, a3, alpha):
                return makeSuperQuadricToroid(x, y, z, a1, a2, a3, alpha, 1.0, 1.0, -math.pi, math.pi, -math.pi, math.pi, 10, 10)
 
x, y, z, rx, ry, rz = 1, 1, 1, 1, 2, 3
cmd.load_cgo(makeEllipsoid(x, y, z, rx, ry, rz), 'ellipsoid-cgo')
x, y, z, rx, ry, rz = 1, 1, 1, 8, 2, 2
cmd.load_cgo(makeToroid(x, y, z, rx, ry, rz, 3), 'toroid-cgo')
