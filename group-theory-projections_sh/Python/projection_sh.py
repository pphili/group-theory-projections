# Import
import numpy as np
import multiprocessing
import functools
import cmath, math
import G8Td as rep #representations of the symmetry group

pi = np.pi

# Functions

# Wigner D Matrices 

def wignerFactorial( s, J, n, m, beta ):
    a1 = int(J + m - s)
    a2 = int(n - m + s) 
    a3 = int(J - n - s)

    if a1 < 0 or a2 < 0 or a3 < 0 :
        return 0
    else:
        
        A1 = math.factorial(a1)
        A2 = math.factorial(a2)
        A3 = math.factorial(a3)

        return ( math.pow(-1, a2)/(A1*math.factorial(s)*A2*A3)*
        math.pow(math.cos(beta/2),int(2*J+m-n-2*s))*
        math.pow(math.sin(beta/2),int(n-m+2*s)) )
        
def d (J, n, m, beta):
    
    maxTerm = int(round(10.0*J))
    
    wdList = [wignerFactorial(x, J, n, m, beta) for x in range(0, maxTerm)]
    
    return math.sqrt(math.factorial(int(J+n))*math.factorial(int(J-n))*math.factorial(int(J+m))*
    math.factorial(int(J-m)))*sum(wdList)
   
def D ( J, n, m, alpha, beta, gamma ):
   return cmath.exp(-1j*n*alpha) * d(J, n, m, beta) * cmath.exp(-1j*m*gamma) 

# Rotation of Spherical Harmonics
def Rot(J, mj, l, alpha, beta, gamma, inv):
    return np.power(inv, l)*np.array([D(J, n-J, mj, alpha, beta, gamma) for n in range(int(2*J+1))])

# Projection
def proj( J, mj, l, i, k ):
    h = 1 * np.trace(rep.G[0])# dimension of the representation
    L = float(len(rep.G)) # dimension of the double group

    s = []
    for jj, A in enumerate(rep.Angles):
        RHO = np.conjugate(rep.G[jj][i,k])
        a1 = A[0]
        a2 = A[1]
        a3 = A[2]
        inv = A[3]

        rot = Rot(J, mj, l, a1, a2, a3, inv)
        s.append(RHO * rot )

    for jj, A in enumerate(rep.Angles):
        RHO = np.conjugate(rep.G[jj+len(rep.Angles)][i,k])
        a1 = A[0] + 2*pi
        a2 = A[1]
        a3 = A[2]
        inv = A[3]

        rot = Rot(J, mj, l, a1, a2, a3, inv)
        s.append(RHO * rot )

    S = sum(s)
    return h / L * S

# Example Valence Band of GaAs
if __name__ == '__main__':

    # Inputs
    # state file
    J = input('J:')
    mj = input('mj:')
    l = input('l:')
    i = input('state to project onto:')
    k = input('state to project from:')

    print('Coefficients printed for states with J={} , l={} with the format :'.format( J, l ))
    print(['mj =' + str(mj) for mj in [-J+n for n in range(int(2*J+1))]])

    print(proj( J, mj, l, i, k ))


   
