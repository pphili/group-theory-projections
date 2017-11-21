# Import
import numpy as np
import multiprocessing
import functools
import reps #represesntations of the symmetry group

# Functions

# Takes file containing state and stores it in appropriate format
def importData(fileName):
    output = []

    for line in open(fileName):
        l = line.split()
        x = [float(elem) for elem in l]
        coord = [round(elem,6) for elem in x[0:3]]
        spin = [x[3]+1j*x[4], x[5]+1j*x[6]]
        s = [coord] + [np.array(spin)]

        output.append(s)

    return output

# Function to perform a rotation on a point
def ROT( pt, r, RS, rho ):

    newr = np.dot(r, np.array(pt[0]))
    newspin = h/L*np.conjugate(rho)*np.dot(RS, pt[1])
    return [[round(elem,6) for elem in newr], newspin]

# Function to get state in an appropriate domain
def cubify( state ):
    pt = state[0][0]
    newstate = []
    for coord in state:
        if pt not in coord[0:3]:
            newstate.append(coord)

    return newstate

# Normalization
def Normalize( state ):
    lc = -float(state[0][0])*2.0
    delta = lc**3/len(state)
    S=sum([np.absolute(x[3][0])**2+np.absolute(x[3][1])**2 for x in state])
    return delta*S

# Projection
def proj( state, cores, names ):
    state.sort()
    # start multiprocessing Pool
    p = multiprocessing.Pool(cores) 

    for ii in range(len(names)):
        projPsi = [np.array([0,0])]*len(state)
        
        # projection
        for jj, rot in enumerate(reps.R):
            rs = reps.wD12[jj]
            RHO = reps.G[jj][ii,ii]

            rotPsi = p.map(functools.partial(ROT, r=rot, RS=rs, rho=RHO), state)
            rotPsi.sort()
            projPsi = [sum(x) for x in zip(projPsi, [y[1] for y in rotPsi])]

        outputPsi = []
        for kk, coord in enumerate(state):
            l = coord[0] + [projPsi[kk]]
            outputPsi.append(l)

        # define wavefunction on a grid for proper normalization
        outputPsi = cubify(outputPsi)
        
        # normalize
        spinors = [x[3] for x in outputPsi]
        N = Normalize(outputPsi)

        finalPsi = []
        for phi in outputPsi:
            nspin = list(np.sqrt(1.0/N)*phi[3])
            x = phi[0:3] + nspin
            finalPsi.append(x)

        # output state
        with open('proj' + str(names[ii]) + '.OUT', 'w') as f:
            for x in finalPsi:
                x = [str(elem) for elem in x]
                f.write('    '.join(x) + '\n')



# Example Valence Band of GaAs
if __name__ == '__main__':

    # Constants
    # Order of the T_d group and \Gamma_8 representation
    h = 4.0 # dimension of the Gamma8 representation
    L = 48.0 # dimension of the T_d double group
    # Integration - GaAs
    # lc = 5.3435 * 2.0 # cubic lattice constant in Bohr radii

    # Inputs
    # state file
    n = input('npts value: ')
    fileName = 'psi3D_npts' + str(n) + '.OUT'

    # Number of cores to be used by the parallelized part of the computation
    nCores = input('number of cores: ')

    # Wavefunction
    Psi = importData(fileName)

    # loop over 4 \Gamma_8 basis states
    labels = ['HHd', 'LHd', 'LHu', 'HHu']
    
    proj(Psi, nCores, labels)
