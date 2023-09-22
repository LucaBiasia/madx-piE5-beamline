import numpy as np
import matplotlib.pyplot as plt
import re

def FieldGrid(fieldmap, skip):
    '''
    This function creates a 6D grid of the given fielmap
    - skip is the number of lines of the header of the fieldmap 
    '''
    fieldmap = fieldmap
    grid = np.loadtxt(fieldmap, dtype=str, skiprows=2, max_rows=1)
    nX = int(*re.findall(r'\d+',grid[4]))
    nY = int(*re.findall(r'\d+',grid[5]))
    nZ = int(*re.findall(r'\d+',grid[6]))
    x, y, z, Bx, By, Bz = np.loadtxt(fieldmap, dtype=float, delimiter='\t', skiprows=skip, usecols=(0,1,2,3,4,5), unpack=True)
    x = np.reshape(x, [nX,nY,nZ])
    y = np.reshape(y, [nX,nY,nZ])
    z = np.reshape(z, [nX,nY,nZ])
    Bx = np.reshape(Bx, [nX,nY,nZ])
    By = np.reshape(By, [nX,nY,nZ])
    Bz = np.reshape(Bz, [nX,nY,nZ])
    return x, y, z, Bx, By, Bz


def Lin(x, a, b, c, d, e, f):
    return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5

def MultipoleFile(filename, X, Z, BY, R0):
    dip_coef, quad_coef, sext_coef, octu_coef, deca_coef, dode_coef = [],[],[],[],[],[]
    coefs = [dip_coef, quad_coef, sext_coef, octu_coef, deca_coef, dode_coef]
    for slice in range(int(len(Z[0,0,:]))):
        x0 = np.flip(np.polyfit(X[:int(R0*100),0,slice]/1000, BY[:int(R0*100),0,slice], len(coefs)-1))
        [coefs[j].append(x0[j]) for j in range(len(coefs))]
    dip_coef = np.append(dip_coef[1:][::-1], dip_coef, axis=0)
    quad_coef = np.append(quad_coef[1:][::-1], quad_coef, axis=0)
    sext_coef = np.append(sext_coef[1:][::-1], sext_coef, axis=0)
    octu_coef = np.append(octu_coef[1:][::-1], octu_coef, axis=0)
    deca_coef = np.append(deca_coef[1:][::-1], deca_coef, axis=0)
    dode_coef = np.append(dode_coef[1:][::-1], dode_coef, axis=0)
    Z = np.append(-Z[0, 0, :][1:][::-1], Z[0, 0, :], axis=0)/1000
    coefs = [Z, dip_coef, quad_coef, sext_coef, octu_coef, deca_coef, dode_coef]
    # Transpose the list of arrays to convert rows to columns
    coefs_transposed = np.transpose(coefs)
    # Save the transposed array to a text file
    np.savetxt(filename, coefs_transposed, fmt='%lf', delimiter='\t')
    return
    

if __name__ == "__main__":
    R0_QSF = 0.2
    R0_HSC = 0.252
    R0_QSB = 0.129
    R0_QSK = 0.125
    FIELDMAP = './fieldmap/feldqsk.dat'
    X, Y, Z, BX, BY, BZ = FieldGrid(FIELDMAP, 7)
    MultipoleFile('QSKmultipoles', X, Z, BY, R0_QSK)