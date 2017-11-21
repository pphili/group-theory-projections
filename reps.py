# Imports
import math
import numpy as np
import cmath
import transforms3d as tr

# This script contains the matrices required to construct and apply the
# group-theoretic projection operators. The matrices stored in 'R' are the
# transfromations for the coordinates, those in 'G' are the desired
# representation of the group and those in 'wD12' are the transformations of the
# spinors. The transformations defining the group are given by 'Angles'. All
# matrices should be stored as numpy arrays. 

pi = np.pi

# Wigner D Matrices 

def wignerFactorial( s, J, n, m, beta ):
    a1 = J + m - s
    a2 = n - m + s 
    a3 = J - n - s

    if a1 < 0 or a2 < 0 or a3 < 0 :
        return 0
    else:
        
        A1 = math.factorial(a1)
        A2 = math.factorial(a2)
        A3 = math.factorial(a3)

        return ( math.pow(-1, a2)/(A1*math.factorial(s)*A2*A3)*
        math.pow(math.cos(beta/2),2*J+m-n-2*s)*
        math.pow(math.sin(beta/2),n-m+2*s) )
        
def d (J, n, m, beta):
    
    maxTerm = int(round(6.0*J))
    
    wdList = [wignerFactorial(x, J, n, m, beta) for x in range(0, maxTerm)]
    
    return math.sqrt(math.factorial(J+n)*math.factorial(J-n)*math.factorial(J+m)*
    math.factorial(J-m))*sum(wdList)
   
def D ( J, n, m, alpha, beta, gamma ):
   return cmath.exp(-1j*n*alpha) * d(J, n, m, beta) * cmath.exp(-1j*m*gamma) 

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

# Construction of rotation matrices for coordinates from the angles
R=[]

for angles in Angles:
    r = angles[3]*tr.euler.euler2mat(angles[0], angles[1], angles[2], axes='rzyz') 
    R.append(np.ndarray.round(r))

#Explicit rotations by 2pi for the double group representations
for angles in Angles:
    r = angles[3]*tr.euler.euler2mat(angles[0]+2*pi, angles[1], angles[2], axes='rzyz') 
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


#Construct matrices for rotations by 2pi for double group
for angles in Angles:
    a = angles[0] + 2*pi
    b = angles[1]
    c = angles[2]
    inv = angles[3]

    # constructing \Gamma_8 representation
    wD = [[math.pow(inv,l)*D(J, n-J, m-J, a, b, c) for m in range(noJStates)] for n in
    range(noJStates)]
    G.append(np.array(wD))

    # constructing matrices describing transformation of spins
    wD = [[D(s, -(n-s), -(m-s), a, b, c) for m in range(noSStates)] for n in
    range(noSStates)]
    wD12.append(np.array(wD))

