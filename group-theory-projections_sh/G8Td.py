# Imports
import math
import numpy as np
import cmath
import transforms3d as tr

# This script contains the matrices required to construct and apply the
# group-theoretic projection operators. The matrices stored in 'G' are the desired
# representation of the group. The transformations defining the group are given 
# by 'Angles'. All matrices should be stored as numpy arrays.  

pi = math.pi

# Angles: first three components are the Euler angles for the rotation and the
# last component is 1 if the symmetry operation contains no inversion and -1 if
# it does

Angles = ( [[0, 0, 0, 1], [pi, pi, 0, 1], [0, pi, 0, 1],
            [pi, 0, 0, 1], [-pi/2,pi/2, pi, 1], [pi/2, pi/2, 0, 1],
            [-pi/2, pi/2, 0, 1], [pi/2, pi/2, pi, 1], [0, pi/2, -pi/2, 1],
            [pi, pi/2, pi/2, 1], [pi, pi/2, -pi/2, 1], [0, pi/2, pi/2, 1],
            [pi/2, pi/2, pi/2, -1], [-pi/2, pi/2, -pi/2, -1],
            [pi/2, pi/2, -pi/2, -1], [-pi/2, pi/2, pi/2, -1],[pi, pi/2, 0, -1],
            [0, pi/2, pi, -1], [pi, pi/2, pi, -1], [0, pi/2, 0, -1], 
            [-pi/2, pi, 0, -1], [pi/2, pi, 0, -1], [-pi/2, 0, 0, -1],
            [pi/2, 0, 0, -1] ] )

#Explicit rotations by 2pi for the double group representations
Angles2 = [[a[0]+2*pi, a[1], a[2], a[3]] for a in Angles]
Angles = Angles + Angles2

# Construction of rotation matrices for coordinates from the angles
R=[]

for angles in Angles:
    r = angles[3]*tr.euler.euler2mat(angles[0], angles[1], angles[2], axes='rzyz') 
    R.append(np.ndarray.round(r))

# Gamma_8 representation of the T_d double group
J = 3./2
l = 1
s = 1./2
noJStates = int(2*J+1)
noSStates = int(2*s+1)

G = []
wD12 = [] 

for angles in Angles:
    a = angles[0]
    b = angles[1]
    c = angles[2]
    inv = angles[3]

    # constructing \Gamma_8 representation
    wD = [[math.pow(inv,l)*D(J, n-J, m-J, a, b, c) for m in
    range(noJStates)] for n in range(noJStates)]
    G.append(np.array(wD))

    # constructing matrices describing transformation of spins
    wD = [[D(s, -(n-s), -(m-s), a, b, c) for m in range(noSStates)] for n in
    range(noSStates)]
    wD12.append(np.array(wD))



