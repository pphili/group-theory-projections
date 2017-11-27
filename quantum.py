import numpy as np

def psiFromFile( fileName, r = 6 ):
    with open(fileName) as f:
        state = [line.split() for line in f]
    state = [[round(float(x), r) for x in l[0:3]] + l[3:5] for l in state]
    state.sort()

    xyz = [[round(float(x), r) for x in l[0:3]] for l in state]
    spinors = [np.array([complex(x) for x in l[3:5]]) for l in state]
    return psi(xyz, spinors)


# Wavefunction class
class psi:
    def __init__( self, coords, wf ):
        self.wf = wf
        self.coords = coords
        
        # computing differential dV
        X = {coord[0] for coord in coords}
        Y = {coord[1] for coord in coords}
        Z = {coord[2] for coord in coords}

        deltaX = (1.+ 1./(len(X)-1))*(max(X) - min(X))
        deltaY = (1.+ 1./(len(Y)-1))*(max(Y) - min(Y))
        deltaZ = (1.+ 1./(len(Z)-1))*(max(Z) - min(Z))

        self.dV = deltaX*deltaY*deltaZ/len(coords)

    # Check if states are defined on the same domain
    def sameDomain( self, state2 ):
        if self.coords == state2.coords:
            return True
        else:
            return False

    # Evaluate matrix elements of the position operator
    def position( self, state2 ):
        x = []
        y = []
        z = []
        for ii, pt in enumerate(self.coords):
            x.append(pt[0]*np.dot(np.conjugate(self.wf[ii]),state2.wf[ii]))
            y.append(pt[1]*np.dot(np.conjugate(self.wf[ii]),state2.wf[ii]))
            z.append(pt[2]*np.dot(np.conjugate(self.wf[ii]),state2.wf[ii]))

        return ( self.dV*sum(x), self.dV*sum(y),self.dV*sum(z) )

    # Evaluate matrix elements of momentum operator
    def p ( self, state2 ):
        X = {coord[0] for coord in self.coords}
        Y = {coord[1] for coord in self.coords}
        Z = {coord[2] for coord in self.coords}

        dX = (1.+ 1./(len(X)-1))*(max(X) - min(X))/len(X)
        dY = (1.+ 1./(len(Y)-1))*(max(Y) - min(Y))/len(Y)
        dZ = (1.+ 1./(len(Z)-1))*(max(Z) - min(Z))/len(Y)

        px = []
        py = []
        pz = []
        
        up = np.array([phi[0] for phi in state2.wf])
        down = np.array([phi[1] for phi in state2.wf])
        
        up = up.reshape((len(X),len(Y),len(Z)))
        down = down.reshape((len(X),len(Y),len(Z)))
        
        momu = np.gradient(up, dX, dY, dZ)
        momd = np.gradient(down, dX, dY, dZ)
        momu = [momentum.reshape(len(self.coords)) for momentum in momu]
        momd = [momentum.reshape(len(self.coords)) for momentum in momd]
        
        PX = [np.array([pu,pd]) for (pu,pd) in zip(momu[0],momd[0])]
        PY = [np.array([pu,pd]) for (pu,pd) in zip(momu[1],momd[1])]
        PZ = [np.array([pu,pd]) for (pu,pd) in zip(momu[2],momd[2])]
        
        for ii, pt in enumerate(self.wf):
            px.append(-1j*np.conjugate(pt).dot(PX[ii]))
            py.append(-1j*np.conjugate(pt).dot(PY[ii]))
            pz.append(-1j*np.conjugate(pt).dot(PZ[ii]))

        return ( self.dV*sum(px), self.dV*sum(py), self.dV*sum(pz) )


    # Matrix elements of the spin operator 
    def spin( self, state2 ):
        paulix = np.array([[0,1],[1,0]])
        pauliy = np.array([[0,-1j],[1j,0]])
        pauliz = np.array([[1,0],[0,-1]])

        x = []
        y = []
        z = []
        for ii, pt in enumerate(self.wf):
            x.append(np.conjugate(pt).dot(paulix).dot(state2.wf[ii]))
            y.append(np.conjugate(pt).dot(pauliy).dot(state2.wf[ii]))
            z.append(np.conjugate(pt).dot(pauliz).dot(state2.wf[ii]))

        return ( self.dV*sum(x), self.dV*sum(y), self.dV*sum(z) )

    # Matrix elements of Orbital Angular Momentum    
    def L ( self, state2 ):
        X = {coord[0] for coord in self.coords}
        Y = {coord[1] for coord in self.coords}
        Z = {coord[2] for coord in self.coords}

        dX = (1.+ 1./(len(X)-1))*(max(X) - min(X))/len(X)
        dY = (1.+ 1./(len(Y)-1))*(max(Y) - min(Y))/len(Y)
        dZ = (1.+ 1./(len(Z)-1))*(max(Z) - min(Z))/len(Y)
        
        Lx = []
        Ly = []
        Lz = []

        up = np.array([phi[0] for phi in state2.wf])
        down = np.array([phi[1] for phi in state2.wf])
        
        up = up.reshape((len(X),len(Y),len(Z)))
        down = down.reshape((len(X),len(Y),len(Z)))

        momu = np.gradient(up, dX, dY, dZ)
        momd = np.gradient(down, dX, dY, dZ)
        momu = [momentum.reshape(len(self.coords)) for momentum in momu]
        momd = [momentum.reshape(len(self.coords)) for momentum in momd]
        
        PX = [np.array([pu,pd]) for (pu,pd) in zip(momu[0],momd[0])]
        PY = [np.array([pu,pd]) for (pu,pd) in zip(momu[1],momd[1])]
        PZ = [np.array([pu,pd]) for (pu,pd) in zip(momu[2],momd[2])]
        
        for ii, pt in enumerate(self.wf):
            Lx.append(-1j*np.conjugate(pt).dot(self.coords[ii][1]*PZ[ii]-self.coords[ii][2]*PY[ii]))
            Ly.append(-1j*np.conjugate(pt).dot(self.coords[ii][2]*PX[ii]-self.coords[ii][0]*PZ[ii]))
            Lz.append(-1j*np.conjugate(pt).dot(self.coords[ii][0]*PY[ii]-self.coords[ii][1]*PX[ii]))


        return ( self.dV*sum(Lx), self.dV*sum(Ly), self.dV*sum(Lz) )

if __name__ == '__main__':
    psi1 = psiFromFile('projHHu.OUT')
    psi2 = psiFromFile('projLHu.OUT')
    # print(psi1.position(psi2))
    # print(psi1.p(psi2))
    print(psi1.L(psi2))
    # print(psi1.spin(psi2))
