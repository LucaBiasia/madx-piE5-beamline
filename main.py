from cpymad.madx import Madx
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import numpy as np
from globals import *
from functions import *
from madx_inputs import *
from my_piE5_beam import *

mx = Madx()


if __name__ == "__main__":
    # Redefine the position along the selected sequence (it is a function of the initial and final positions and the direction of the propagation)
    for key in position.keys():
        position[key] = (position[key]-POS_IN)*(1+REVERSE)/2 + (POS_OUT-position[key])*(1-REVERSE)/2
    
    slicing, dip_spline, quad_spline, sextu_spline, octu_spline, deca_spline, dodeca_spline, solenoids_spline = CreateSplines(REVERSE, POS_IN, POS_OUT, DX, P0, FIELDMAP_PATH)
    # Defining variables
    mx.input(f'''
    sigmax := {SIGMAX};     // m
    sigmaxp := {SIGMAXP};     // rad
    sigmay := {SIGMAY};     // m
    sigmayp := {SIGMAYP};     // rad
    rhox = {RHOX};    // x correlation
    rhoy = {RHOY};    // y correlation
    ex := {ex};
    ey := {ey};
    alfax = {alfax};
    alfay = {alfay};

    // BEAM DEFINITION----------------------------------------------------------
    beam, particle=POSMUON, pc={P0}, EX={ex}, EY={ey}, bunched=False;
    ''')
    # Defining sequence
    mx.input(f'''
    // ELEMENT DEFINITION-------------------------------------------------------
    AHSW : matrix, {RM_AHSW41}, {TM_AHSW41};

    ASC : sbend, angle={np.radians(75)}, l={L_ASC}, e1={-np.radians(7.5)}, e2={-np.radians(7.5)}, hgap=0.3, thick=true;

    AST : sbend, angle={-np.radians(47.5)}, l={L_AST}, e1={np.radians(19)}, h1=-0.3, h2=0, hgap=0.12, fint=0.355, thick=true;

    MP : quadrupole,    l={DX}, thick=true;

    SOL : solenoid,     l={DX}, thick=true;
    // SEQUENCE DEFINITION------------------------------------------------------
    piE5: sequence,  l={L_TOT}, refer=center; // with refer=entry positions are given from the entry of the element (otherwise =exit or =center)
    ''')

    exit1 = 0
    exit2 = 0
    # exit3 = 0
    for i in range(len(slicing)):
        if ((slicing[i]+DX/2 > position['poszASC41']-L_ASC/2-DX) and (slicing[i]-DX/2 < position['poszASC41']+L_ASC/2+DX) and (exit1==0)):
            mx.input('''
                ASC41 : ASC, at=%f;
            '''%(position['poszASC41']))
            exit1 += 1
        elif ((slicing[i]+DX/2 > position['poszAST41']-L_AST/2) and (slicing[i]-DX/2 < position['poszAST41']+L_AST/2) and (exit2==0)):
            mx.input('''
                AST41 : AST, at=%f;
            '''%(position['poszAST41']))
            exit2 += 1
        elif (not (((slicing[i]+DX/2 > position['poszAST41']-L_AST/2) and (slicing[i]-DX/2 < position['poszAST41']+L_AST/2)) or ((slicing[i]+DX/2 > position['poszASC41']-L_ASC/2-DX) and (slicing[i]-DX/2 < position['poszASC41']+L_ASC/2+DX)) or (slicing[i] > (position['posQSF41']+0.5)))):
            mx.input('''
                SOL%i : SOL, at=%f, ks=%f, knl={%f, %f, %f, %f, %f, %f};
            ''' %(slicing[i]*1000, slicing[i], solenoids_spline(slicing[i]), dip_spline(slicing[i])*DX, quad_spline(slicing[i])*DX, sextu_spline(slicing[i])*DX, octu_spline(slicing[i])*DX, deca_spline(slicing[i])*DX, dodeca_spline(slicing[i])*DX))
        # elif ((slicing[i] > (position['posQSF41']+0.5)) and (exit3==0)):
        #     mx.input('''
        #     AHSW41 : AHSW, at=%f;
        #     '''%(position['posfrontarcAHSW41']))
        #     exit3 += 1
    
    mx.input(end_sequence)

    ##########################################  PTC TWISS ##############################################
    mx.input(ptc_twiss)
    ptc_twiss = mx.table.ptc_twiss.dframe()
    # Twiss parameters
    moments = mx.table['moments']
    betx = moments['mu200000']
    meanx2TWISS = moments['mu200000']
    bety = moments['mu002000']
    rhoxTWISS = moments['mu110000']
    rhoyTWISS = moments['mu001100']
    sigmaxpTWISS = moments['mu020000']
    sigmaypTWISS = moments['mu000200']
    meanxTWISS = moments['mu100000']
    meanyTWISS = moments['mu001000']



    ########################################## PTC TRACKING #############################################
    # MULTIVARIATE SAMPLING
    x_ps, px_ps, y_ps, py_ps = BeamSampling(N_SAMPLE, N_G4BL, SIGMAX, SIGMAXP, SIGMAY, SIGMAYP, RHOX, RHOY)

    N_SAMPLE = len(x_ps)
    np.random.seed(0)
    ptot = np.random.normal(P0, DP0, N_SAMPLE) #Gev/c

    # Evaluate energy error
    MUON_MASS = 0.1056583755
    pt_ps = (ptot-P0)/P0

    # Initializing the tracking procedure
    mx.input(ptc_create)

    # Setting the sampled particles
    for i in range(N_SAMPLE):
        mx.input(f'''
        ptc_start, x={x_ps[i]}, px={px_ps[i]}, y={y_ps[i]}, py={py_ps[i]}, t=0, pt={pt_ps[i]}; // pt is the energy difference divided by the reference momentum
        ''')

    # Setting the observation points
    for ele in mx.elements.keys():
        mx.input(f'''ptc_observe, place={ele}''')

    # Tracking the particles through the sequence
    mx.input(ptc_track)

    ptc_output = mx.table.trackone.dframe()
    ptc_lost = mx.table.trackloss.dframe()

    x = ptc_output[abs(ptc_output.s -L_TOT) < 0.001].x*1000
    y = ptc_output[abs(ptc_output.s -L_TOT) < 0.001].y*1000
    px = ptc_output[abs(ptc_output.s -L_TOT) < 0.001].px
    py = ptc_output[abs(ptc_output.s -L_TOT) < 0.001].py
    pt = ptc_output[abs(ptc_output.s -L_TOT) < 0.001].pt
    ptot_f = np.sqrt((pt*P0+np.sqrt(P0**2+MUON_MASS**2))**2-MUON_MASS**2)
    pz = ptot_f/np.sqrt(1+px**2+py**2)*1000

    #-------------------------------------------------- G4BL
    # CREATING A ROOT NTUPLE FOR G4BEAMLINE WITH THE SAMPLED TRACKS
    # PATH = '/home/lucabiasia/Scrivania/piE5/'
    # CreateNTuple(x_ps, px_ps, y_ps, py_ps, ptot, PATH=PATH)

    # RUNNING G4BL
    # profile_name = 'TESTMADX_profile.txt'
    # command = 'g4bl /home/lucabiasia/Scrivania/piE5/TESTMADX.g4bl ' + f"QSK41cur={QSK41['current']} " + f"QSK42cur={QSK42['current']} " + f"QSK43cur={QSK43['current']} " +  f"QSKSCALE_={-abs(QSK41['scaling'])} " + f'FIELDMAPQSK={QSK41["fieldmap"]} ' + f'FIELDMAPBTS={BTS["fieldmap"]} ' + f'FIELDMAPCOBRA={COBRA["fieldmap"]} '
    # command += f"profilename=/home/lucabiasia/Scrivania/piE5/profiles/{profile_name}"
    # os.system(command)

    # zloop, mux, sigx, muy, sigy = np.loadtxt(f'/home/lucabiasia/Scrivania/single_element/profiles/{profile_name}', delimiter=' ', usecols=(0, 2, 3, 4, 5), unpack=True)


    # g4bl_output = uproot.open(PATH+'g4beamline.root')
    # virtual_detector = g4bl_output['/VirtualDetector/VirtualDetector1']
    # x_g4 = virtual_detector['x'].array(library="np")
    # y_g4 = virtual_detector['y'].array(library="np")
    # Px_g4 = virtual_detector['Px'].array(library="np")
    # Py_g4 = virtual_detector['Py'].array(library="np")
    # Pz_g4 = virtual_detector['Pz'].array(library="np")

    # xp_g4 = Px_g4/Pz_g4 # rad
    # yp_g4 = Py_g4/Pz_g4 # rad



    #----------------------------------------- PHACE SPACE AND BEAM SPOT COMPARISON BETWEEN MAD-X AND G4BL -------------------------------------------

    plt.figure('Particle tracks', figsize=(10, 6))
    plt.title(f'Particle tracks in piE5 from {np.round(POS_IN,2)}m to {np.round(POS_OUT)}: N = {N_SAMPLE}', fontsize=16)
    for i in range(int(N_SAMPLE*2)+1):
        i=i+1
        # Tracks
        plt.plot(ptc_output[ptc_output.number == i].s, np.abs(ptc_output[ptc_output.number == i].x)*1000, color="b", alpha=0.1)
        plt.plot(ptc_output[ptc_output.number == i].s, -np.abs(ptc_output[ptc_output.number == i].y)*1000, color="r", alpha=0.1)
    # MAD-X Twiss observation points
    plt.plot(ptc_twiss['s'], (np.sqrt(betx))*1000,'o-b',label='$\\sigma_x$ MADX')
    plt.plot(ptc_twiss['s'], -(np.sqrt(bety))*1000,'o-r',label='$\\sigma_x$ MADX')

    # mean of the envelope
    # plt.plot(ptc_twiss['s'], meanxTWISS*1000)
    # plt.plot(ptc_twiss['s'], meanyTWISS*1000)
    
    # G4BL envelope
    # plt.plot(zloop/1000, sigx, '-', linewidth=3,label='$\\sigma_x$ G4BL')
    # plt.plot(zloop/1000, -sigy, '-', linewidth=3,label='$\\sigma_y$ G4BL')
    # plt.plot(zloop/1000, mux, '-', linewidth=3,label='$\\mu_x$ G4BL')
    # plt.plot(zloop/1000, -muy, '-', linewidth=3,label='$\\mu_y$ G4BL')

    for ele in Sequence:
        ele['position'] = position[f'pos{ele["name"]}']
        left, bottom, width, height = (ele['position']-ele['L']/2, -ele['R0']*1000, ele['L'], ele['R0']*2*1000)
        if (ele['type'] == 'QSF'):
            rect=mpatches.Rectangle((left,bottom),width,height,
                                    alpha=0.1,
                                    facecolor="red")
        elif (ele['type'] == 'QSB'):
            rect=mpatches.Rectangle((left,bottom),width,height,
                                    alpha=0.1,
                                    facecolor="red")
        elif (ele['type'] == 'QSK'):
            rect=mpatches.Rectangle((left,bottom),width,height,
                                    alpha=0.1,
                                    facecolor="red")
        else:
            rect=mpatches.Rectangle((left,bottom),width,height,
                                    alpha=0.1,
                                    facecolor="green")
        plt.gca().add_patch(rect)

    left, bottom, width, height = (position['poszAST41']-L_AST/2, -R0_AST*1000, L_AST, R0_AST*2*1000)
    rect=mpatches.Rectangle((left,bottom),width,height,
                        alpha=0.1,
                    facecolor="blue")
    plt.gca().add_patch(rect)

    left, bottom, width, height = (position['poszASC41']-L_ASC/2, -R0_ASC*1000, L_ASC, R0_ASC*2*1000)
    rect=mpatches.Rectangle((left,bottom),width,height,
                        alpha=0.1,
                    facecolor="blue")
    plt.gca().add_patch(rect)

    left, bottom, width, height = (position['posBTS']-L_BTS/2, -R0_BTS*1000, L_BTS, R0_BTS*2*1000)
    rect=mpatches.Rectangle((left,bottom),width,height,
                        alpha=0.1,
                    facecolor="grey")
    plt.gca().add_patch(rect)

    left, bottom, width, height = (position['posCOBRA']-L_COBRA/2, -R0_COBRA*1000, L_COBRA, R0_COBRA*2*1000)
    rect=mpatches.Rectangle((left,bottom),width,height,
                        alpha=0.1,
                    facecolor="green")
    plt.gca().add_patch(rect)

    # plt.text(pos41-L_QSK/2, -150, 'QSK41', rotation='horizontal', fontsize=9)
    # plt.text(pos42-L_QSK/2, -150, 'QSK42', rotation='horizontal', fontsize=9)
    # plt.text(pos43-L_QSK/2, -150, 'QSK43', rotation='horizontal', fontsize=9)
    # plt.text(position['posBTS'], -150, 'BTS', rotation='horizontal', fontsize=14)
    # plt.text(position['posCOBRA'], -150, 'COBRA', rotation='horizontal', fontsize=14)
    plt.vlines(position['posMEGCOL'], -max(np.sqrt(betx)*1000)*1.2, max(np.sqrt(betx)*1000)*1.2, linestyles='dashed')
    plt.xlabel('Position along the sequence (m)')
    plt.ylabel('Distance from axis (mm)')
    plt.legend()
    plt.grid(alpha=0.2)

    # print(np.sqrt(betx[-1]))
    # print(rhoxTWISS[-1]/(np.sqrt(betx[-1])*np.sqrt(sigmaxpTWISS[-1])))
    # print(np.sqrt(sigmaxpTWISS[-1]))
    # print(np.sqrt(bety[-1]))
    # print(rhoyTWISS[-1]/(np.sqrt(bety[-1])*np.sqrt(sigmaypTWISS[-1])))
    # print(np.sqrt(sigmaypTWISS[-1]))
    # print(ptc_twiss.keys())
    # sx = np.sqrt(ptc_twiss['betx']*ex)
    # sxp = np.sqrt((1 + ptc_twiss['alfx']**2)/ptc_twiss['betx']*ex)
    # print(-(ptc_twiss['alfx']*ex/sx/sxp)[-1])




    plt.show()


# NORMALIZED PHASE SPACE
    # twiss = mx.table.twiss.dframe() # CREATING DATAFRAME WITH TWISS CALCULATION
    # Rx = (1/np.sqrt(twiss['betx'][-1]))*np.array(((1, 0), (twiss['alfx'][-1], twiss['betx'][-1])))
    # Ry = (1/np.sqrt(twiss['bety'][-1]))*np.array(((1, 0), (twiss['alfy'][-1], twiss['bety'][-1])))
    # xx, pxx = np.dot(Rx, [ptc_output[abs(ptc_output.s -L_TOT) < 0.001].x*1000, ptc_output[abs(ptc_output.s -L_TOT) < 0.001].px*1000])
    # yy, pyy = np.dot(Ry, [ptc_output[abs(ptc_output.s -L_TOT) < 0.001].y*1000, ptc_output[abs(ptc_output.s -L_TOT) < 0.001].py*1000])

    # plt.figure(figsize=(12, 10))
    # plt.suptitle(f'from QSK41 to COBRA center', fontsize=16)
    # plt.subplot(2,2,1)
    # plt.title('Normalized X phace space')
    # plt.plot(x, px, '.', color="b", alpha=0.3, label='MAD-X')
    # plt.xlabel('position (mm)')
    # plt.ylabel('divergence (mrad)')
    # plt.grid(alpha=0.3)
    # plt.legend()

    # plt.subplot(2,2,3)
    # plt.title('Normalized Y phace space')
    # plt.plot(y, py, '.', color="b", alpha=0.3, label='MAD-X')
    # plt.xlabel('position (mm)')
    # plt.ylabel('divergence (mrad)')
    # plt.grid(alpha=0.3)
    # plt.legend()

    # plt.subplot(1,2,2)
    # plt.title('Beam spot')
    # plt.plot(ptc_output[abs(ptc_output.s -L_TOT) < 0.001].x*1000, ptc_output[abs(ptc_output.s -L_TOT) < 0.001].y*1000, '.', color="b", alpha=0.3, label='MAD-X')
    # plt.xlabel('x (mm)')
    # plt.ylabel('y (mm)')
    # plt.grid(alpha=0.3)
    # plt.legend()


