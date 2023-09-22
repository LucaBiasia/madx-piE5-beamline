from cpymad.madx import Madx
import numpy as np
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import uproot
import os
import re
from globals import *

######################################### FUNCTIONS ##############################################

def SampleEllipse(xmax, ymax, phi, theta, x0=0, y0=0):
    '''
    Function that extrapolates the minor and major radius of a generic tilted ellipse (e.g. Phase space) knowing only the maximum of x and y.
    Parameters
    ----------
        xmax : maximum of the x axis
        ymax : maximum of the y axis
        phi : angle of the tilted ellipse = arctan(2*alpha/(gamma-beta))/2
        theta : rotating angle around the center (from 0 to 2*pi)
        x0 : x coordinate of the centre of the ellipse
        y0 : y coordinate of the centre of the ellipse
    Returns
    ---------
        Coordinate x and y of a generic point of the ellipse.
    '''
    Ra = np.sqrt((xmax**2/np.cos(phi)**2-ymax**2*np.sin(phi)**2/np.cos(phi)**4)/(1-np.tan(phi)**4))
    Rb = np.sqrt((ymax**2-Ra**2*np.sin(phi)**2)/np.cos(phi)**2)
    x = Ra*np.cos(theta)*np.cos(phi) - Rb*np.sin(theta)*np.sin(phi) + x0
    y = Ra*np.cos(theta)*np.sin(phi) + Rb*np.sin(theta)*np.cos(phi) + y0
    return x, y


def CreateNTuple(x, xp, y, yp, Ptot, PATH):
    '''
    This function creates a Root NTuple already in a G4beamline format, given the sampled particles variables.
    Parameters
    ----------
        x : sigma_x of the particles (m)
        xp : x divergence of the particles (rad)
        y : sigma_y of the particles (m)
        yp : y divergence of the particles (rad)
        Ptot : Total momentum of the sampled particles (Gev/c)
        PATH : path to the folder to create the ROOT NTuple
    Returns
    ---------
        A Root NTuple called InitPhaseSpace.root in the same folder
    Also z and t should be included (for a non-continuous beam)
    
    N.B. This function requires TreeToNTuple.cpp file
    '''
    N = len(x)

    x = x*1000  #mm
    y = y*1000  #mm
    Pz = Ptot*1000/np.sqrt(1+xp**2+yp**2) #MeV/c
    Px = Pz*xp #MeV/c
    Py = Pz*yp #MeV/c

    z = np.full(N, 0)
    t = np.full(N, 0)
    PDGid = np.full(N, -13)
    EventID = np.arange(1, N+1)
    TrackID = np.full(N, 0)
    ParentID = np.full(N, 0)
    Weight = np.full(N, 1)

    dict = {
        'x' : x.astype('float32'),
        'y' : y.astype('float32'),
        'z' : z.astype('float32'),
        'Px' : Px.astype('float32'),
        'Py' : Py.astype('float32'),
        'Pz' : Pz.astype('float32'),
        't' : t.astype('float32'),
        'PDGid' : PDGid.astype('float32'),
        'EventID' : EventID.astype('float32'),
        'TrackID' : TrackID.astype('float32'),
        'ParentID' : ParentID.astype('float32'),
        'Weight' : Weight.astype('float32')
    }

    tree = uproot.recreate("Init_tmp.root")
    tree['Init_tmp'] = dict
    tree.close()
    command = 'root -q ' + PATH + 'TreeToNTuple.cpp'
    os.system(command)
    os.remove(PATH + "Init_tmp.root")
    return 0

def FieldGrid(fieldmap, skip):
    '''
    This function creates a 6D grid of the given fieldmap
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

def FieldGridAHSW(fieldmap, skip):
    '''
    This function creates a 6D grid of the given fieldmap
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


def BTSGrid(fieldmap, skip):
    '''
    This function creates a 6D grid of the given fieldmap
    - skip is the number of lines of the header of the fieldmap 
    '''
    fieldmap = fieldmap
    grid = np.loadtxt(fieldmap, dtype=str, skiprows=2, max_rows=1)
    nR = int(*re.findall(r'\d+',grid[2]))
    nZ = int(*re.findall(r'\d+',grid[3]))
    r, z, Br, Bz = np.loadtxt(fieldmap, dtype=float, delimiter='\t', skiprows=skip, usecols=(0,1,2,3), unpack=True)
    r = np.reshape(r, [nZ,nR])
    z = np.reshape(z, [nZ,nR])
    Br = np.reshape(Br, [nZ,nR])
    Bz = np.reshape(Bz, [nZ,nR])
    return r, z, Br, Bz

def EffectiveLength(fieldmap, skip):
    '''
    This function computes the effective length of a magnet for a given fieldmap
    - skip is the number of lines of the header of the fieldmap 
    '''
    _, _, Z, _, BY, _ = FieldGrid(fieldmap, skip)
    z = Z[1,0,:]
    By = BY[1,0,:]
    By_max = np.max(np.abs(By))
    f = interp1d(z, np.abs(By), kind='cubic')
    I = integrate.quad(f, z[0], z[-1])
    L = I[0]*2/By_max
    return np.round(L, 2)

def InverseDrift(x_f, xp_f, L):
    '''
    This function computes a backwards drift of length L, given the final phase space
    '''
    Sigma = np.array([[x_f**2, x_f*xp_f],
            [xp_f*x_f, xp_f**2]])
    M_drift = np.array([[1, L],
            [0, 1]])
    M_drift_inv = np.linalg.inv(M_drift)
    Sigma_in = np.dot(np.dot(Sigma, M_drift.T), M_drift_inv) #np.linalg.multi_dot(M_drift_inv, Sigma, M_drift.T)
    x_in = Sigma_in[0,0]
    xp_in = Sigma_in[1,1]
    return x_in, xp_in

def SumFunc(*args):
    '''
    This function computes the summation of an arbitrary number of functions (e.g. splines)
    '''
    def compute(x):
        func = 0
        for f in args:
            func += f(x)
        return func
    return compute

def CreateSplines(REVERSE, POS_IN, POS_OUT, DX, P0, FIELDMAP_PATH):
    L_TOT = POS_OUT - POS_IN
    f = 0.3/P0    # conversion factor from T/m to m^-2 (1/(B*rho))
    # quadrupoles and sextupoles
    slicing = np.arange(DX/2, L_TOT-DX, DX)
    dip_spline = []
    quad_spline = []
    sextu_spline = []
    octu_spline = []
    deca_spline = []
    dodeca_spline = []
    for element in Sequence:
        z, dip, quad, sextu, octu, deca, dodeca = np.loadtxt('multipoles/'+element['type']+'multipoles', delimiter='\t', unpack=True)
        # Calculating the strengh of the multipole coefficients
        dip_strength = dip*f*element['current']*element['scaling']*REVERSE
        quad_strength = quad*f*element['current']*element['scaling']
        sextu_strength = sextu*f*element['current']*element['scaling']*2*REVERSE    # sextupoles are not symmetrical under x -> -x
        octu_strength = octu*f*element['current']*element['scaling']*3
        deca_strength = deca*f*element['current']*element['scaling']*4*REVERSE
        dodeca_strength = dodeca*f*element['current']*element['scaling']*5
        # Creating multipole splines for quadrupoles and sextupoles
        z = z + element['position']
        dip_spline.append(interp1d(z, dip_strength, kind='cubic', bounds_error=False, fill_value=0))
        quad_spline.append(interp1d(z, quad_strength, kind='cubic', bounds_error=False, fill_value=0))
        sextu_spline.append(interp1d(z, sextu_strength, kind='cubic', bounds_error=False, fill_value=0))
        octu_spline.append(interp1d(z, octu_strength, kind='cubic', bounds_error=False, fill_value=0))
        deca_spline.append(interp1d(z, deca_strength, kind='cubic', bounds_error=False, fill_value=0))
        dodeca_spline.append(interp1d(z, dodeca_strength, kind='cubic', bounds_error=False, fill_value=0))
    dip_spline = SumFunc(*dip_spline)
    quad_spline = SumFunc(*quad_spline)
    sextu_spline = SumFunc(*sextu_spline)
    octu_spline = SumFunc(*octu_spline)
    deca_spline = SumFunc(*deca_spline)
    dodeca_spline = SumFunc(*dodeca_spline)

    # Creating splines for BTS
    BTSmap = FIELDMAP_PATH + '/fieldBTS.dat'
    _, Z, _, BZ = BTSGrid(BTSmap, skip=5)
    zBTS = Z[:,0]/1000 + position['posBTS']
    zBTS = np.append(-zBTS[1:][::-1], zBTS, axis=0)
    BzBTS = BZ[:,0] * f * BTS['scaling']
    BzBTS = np.append(BzBTS[1:][::-1], BzBTS, axis=0)
    BTS_spl = interp1d(zBTS, BzBTS, kind='cubic', bounds_error=False, fill_value=0)

    # Creating splines for COBRA
    COBRAmap = FIELDMAP_PATH + '/COBRAmap.dat'
    _, _, Z, _, _, BZ = FieldGrid(COBRAmap, skip=4)
    zCOBRA = Z[45,45,:]/1000 + position['posCOBRA']
    BzCOBRA = BZ[45,45,:] * f * COBRA['scaling']
    COBRA_spl = interp1d(zCOBRA, BzCOBRA, kind='cubic', bounds_error=False, fill_value=0)
    solenoids_spline = SumFunc(BTS_spl, COBRA_spl)

    return slicing, dip_spline, quad_spline, sextu_spline, octu_spline, deca_spline, dodeca_spline, solenoids_spline

def BeamSampling(N_SAMPLE, N_G4BL, SIGMAX, SIGMAXP, SIGMAY, SIGMAYP, RHOX, RHOY):
    ex = SIGMAX*SIGMAXP*np.sqrt(1-RHOX**2)      # x emittance
    ey = SIGMAY*SIGMAYP*np.sqrt(1-RHOY**2)      # y emittance
    alfax = -RHOX*SIGMAX*SIGMAXP/ex
    alfay = -RHOY*SIGMAY*SIGMAYP/ey
    betax = SIGMAX**2/ex
    betay = SIGMAY**2/ey
    gammax = (1 + alfax**2)/betax
    gammay = (1 + alfay**2)/betay
    mean = [0, 0, 0, 0]
    cov = np.array([[betax*ex, -alfax*ex, 0, 0],
        [-alfax*ex, gammax*ex, 0, 0],
        [0, 0, betay*ey, -alfay*ey],
        [0, 0, -alfay*ey, gammay*ey]])
    np.random.seed(1)
    phace_space4D = np.random.multivariate_normal(mean, cov, size=N_SAMPLE)
    phace_space4D_ = phace_space4D[:N_G4BL]
    x_ps = phace_space4D_[:,0]
    px_ps = phace_space4D_[:,1]
    y_ps = phace_space4D_[:,2]
    py_ps = phace_space4D_[:,3]

    N_SAMPLE = len(x_ps)
    np.random.seed(0)

    temp = []
    cov = np.linalg.inv(cov)
    for p in phace_space4D:
        temp.append(np.dot(p.transpose(), cov))
        temp[-1] = np.dot(temp[-1].transpose(), p)
    temp = np.array(temp)
    phace_space4D = phace_space4D[(temp < 1)*(temp>0.99)]
    x_ps = phace_space4D[:,0]
    px_ps = phace_space4D[:,1]
    y_ps = phace_space4D[:,2]
    py_ps = phace_space4D[:,3]
    return x_ps, px_ps, y_ps, py_ps
