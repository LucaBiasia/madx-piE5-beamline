#######################################      GLOBALS      ##########################################
import numpy as np


PATH = '/home/lucabiasia/Scrivania/piE5/piE5_for_TRIUMF_noBeam/'

##################################### ELEMENTS POSITIONS ########################################
pos_file = np.loadtxt(PATH+'Positions.txt', skiprows=1, dtype=str)[:,2]

position = {}
for line in pos_file:
    name, val = line.split('=')
    position[name] = float(val)/1000 #position in meters



##################################### ELEMENTS CURRENTS ##########################################
curr_file = np.loadtxt(PATH+'CurrentsCMBL2018.txt', skiprows=1, dtype=str)[:,2]

current = {}
for line in curr_file:
    name, val = line.split('=')
    current[name] = float(val) #currents are already in Tesla



##################################### SCALING FACTORS ###########################################
scale_file = np.loadtxt(PATH+'Scalefactors.txt', skiprows=1, max_rows=21, dtype=str)[:,2]

scaling = {}
for line in scale_file:
    name, val = line.split('=')
    scaling[name] = float(eval(val))


######################################## MAGNET PARAMETERS #########################################
# Scaling factors
SCALE_QSF = scaling['scaleQSF']
SCALE_HSC = scaling['scaleHSC']*-1
SCALE_QSB = scaling['scaleQSB']
SCALE_QSK = scaling['scaleQSK']*-1
SCALE_COBRA = 1.
SCALE_BTS = -198.6/200

# Conversion factors of piE5 elements (T/A)
CF_QSF = 9.0516*1e-4
CF_HSC = 9.1580*1e-4
CF_AST = 12.727*1e-4
CF_ASC = 15.356*1e-4
CF_QSB = 17.5878*1e-4
CF_QSK = 18.0255*1e-4

# Radii of magnets
R0_QSF = 0.2
R0_HSC = 0.252
R0_ASC = 0.3
R0_AST = 0.12
R0_QSB = 0.129
R0_QSK = 0.125
R0_BTS = 0.25
R0_COBRA = 0.32

# Effective lengths of magnets (m)
L_QSF = 0.382
L_HSC = 0.326
L_ASC = 0.640
L_AST = 0.926
L_QSB = 0.314
L_QSK = 0.390
L_BTS = 2.808
L_COBRA = 2.840

# Currents of magnets (A)
QSF41_CURRENT = -93.61     #-93.251
HSC41_CURRENT = 2.18       #-19.06
QSF42_CURRENT = 93.92      #82.6211
QSF43_CURRENT = -78.25     #-64.4526
HSC42_CURRENT = -36.72     #-29.8842
QSF44_CURRENT = 75.93      #81.8146
QSF45_CURRENT = 63.31      #84.9998
HSC43_CURRENT = -13.79     #44.9989
QSF46_CURRENT = -93.03     #-82.662
QSF47_CURRENT = 105.52     #95.1497
HSC44_CURRENT = -7.37      #-10
QSF48_CURRENT = -60.47     #-50.5661
QSB41_CURRENT = 49.55      #37.76
QSB42_CURRENT = -60.32     #-48
QSB43_CURRENT = 34.95      #15
QSK41_CURRENT = 31.0       #22.17
QSK42_CURRENT = -49.54     #-41.16
QSK43_CURRENT = 45.23      #33.04


FIELDMAP_PATH = '/home/lucabiasia/Scrivania/fieldmaps/'
# Straight line
QSF41 = {"name" : "QSF41", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF41'], 'R0' : R0_QSF, 'current' : QSF41_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
HSC41 = {"name" : "HSC41", "type" : "HSC", 'fieldmap' : FIELDMAP_PATH + 'feldhsc.dat', 'position': position['posHSC41'], 'R0' : R0_HSC, 'current' : HSC41_CURRENT, 'scaling' : SCALE_HSC, 'L' : L_HSC}
QSF42 = {"name" : "QSF42", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF42'], 'R0' : R0_QSF, 'current' : QSF42_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
QSF43 = {"name" : "QSF43", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF43'], 'R0' : R0_QSF, 'current' : QSF43_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
HSC42 = {"name" : "HSC42", "type" : "HSC", 'fieldmap' : FIELDMAP_PATH + 'feldhsc.dat', 'position': position['posHSC42'], 'R0' : R0_HSC, 'current' : HSC42_CURRENT, 'scaling' : SCALE_HSC, 'L' : L_HSC}
QSF44 = {"name" : "QSF44", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF44'], 'R0' : R0_QSF, 'current' : QSF44_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
QSF45 = {"name" : "QSF45", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF45'], 'R0' : R0_QSF, 'current' : QSF45_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
HSC43 = {"name" : "HSC43", "type" : "HSC", 'fieldmap' : FIELDMAP_PATH + 'feldhsc.dat', 'position': position['posHSC43'], 'R0' : R0_HSC, 'current' : HSC43_CURRENT, 'scaling' : SCALE_HSC, 'L' : L_HSC}
QSF46 = {"name" : "QSF46", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF46'], 'R0' : R0_QSF, 'current' : QSF46_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
QSF47 = {"name" : "QSF47", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF47'], 'R0' : R0_QSF, 'current' : QSF47_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
HSC44 = {"name" : "HSC44", "type" : "HSC", 'fieldmap' : FIELDMAP_PATH + 'feldhsc.dat', 'position': position['posHSC44'], 'R0' : R0_HSC, 'current' : HSC44_CURRENT, 'scaling' : SCALE_HSC, 'L' : L_HSC}
QSF48 = {"name" : "QSF48", "type" : "QSF", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSF48'], 'R0' : R0_QSF, 'current' : QSF48_CURRENT, 'scaling' : SCALE_QSF, 'L' : L_QSF}
# Triplet-I
QSB41 = {"name" : "QSB41", "type" : "QSB", 'fieldmap' : FIELDMAP_PATH + 'feldqsb.dat', 'position': position['posQSB41'], 'R0' : R0_QSB, 'current' : QSB41_CURRENT, 'scaling' : SCALE_QSB, 'L' : L_QSB}
QSB42 = {"name" : "QSB42", "type" : "QSB", 'fieldmap' : FIELDMAP_PATH + 'feldqsb.dat', 'position': position['posQSB42'], 'R0' : R0_QSB, 'current' : QSB42_CURRENT, 'scaling' : SCALE_QSB, 'L' : L_QSB}
QSB43 = {"name" : "QSB43", "type" : "QSB", 'fieldmap' : FIELDMAP_PATH + 'feldqsb.dat', 'position': position['posQSB43'], 'R0' : R0_QSB, 'current' : QSB43_CURRENT, 'scaling' : SCALE_QSB, 'L' : L_QSB}
# Triplet-II
QSK41 = {"name" : "QSK41", "type" : "QSK", 'fieldmap' : FIELDMAP_PATH + 'feldqsk.dat', 'position': position['posQSK41'], 'R0' : R0_QSK, 'current' : QSK41_CURRENT, 'scaling' : SCALE_QSK, 'L' : L_QSK}
QSK42 = {"name" : "QSK42", "type" : "QSK", 'fieldmap' : FIELDMAP_PATH + 'feldqsk.dat', 'position': position['posQSK42'], 'R0' : R0_QSK, 'current' : QSK42_CURRENT, 'scaling' : SCALE_QSK, 'L' : L_QSK}
QSK43 = {"name" : "QSK43", "type" : "QSK", 'fieldmap' : FIELDMAP_PATH + 'feldqsf.dat', 'position': position['posQSK43'], 'R0' : R0_QSK, 'current' : QSK43_CURRENT, 'scaling' : SCALE_QSK, 'L' : L_QSK}
# BTS and COBRA
BTS = {"name" : "BTS", "type" : "BTS", 'fieldmap' : FIELDMAP_PATH + 'fieldBTS.dat', 'position': position['posBTS'], 'R0' : R0_BTS, 'current' : 1, 'scaling' : SCALE_BTS, 'L' : L_BTS}
COBRA = {"name" : "COBRA", "type" : "COBRA", 'fieldmap' : FIELDMAP_PATH + 'COBRAmap.dat', 'position': position['posCOBRA'], 'R0' : R0_COBRA, 'current' : 1, 'scaling' : SCALE_COBRA, 'L' : L_COBRA}

Sequence = [QSF41, HSC41, QSF42, QSF43, HSC42, QSF44, QSF45, HSC43, QSF46, QSF47, HSC44, QSF48, QSB41, QSB42, QSB43, QSK41, QSK42, QSK43] #, BTS, COBRA

# Rmatrix and Tmatrix of AHSW41
RM_AHSW41 = ''' rm11=-0.60096, rm12=0.03516, rm13=0, rm14=0, rm15=0, rm16=49321,
    rm21=-18.52476, rm22=-0.57297, rm23=0, rm24=0, rm25=0, rm26=5.26292,
    rm31=0, rm32=0, rm33=2.71348, rm34=0.26018, rm35=0, rm36=0,
    rm41=0, rm42=0, rm43=24.45407, rm44=2.68069, rm45=0, rm46=0
    '''

TM_AHSW41 = ''' tm111=0.002966,      tm112=0.00079260,    tm113=0,          tm114=0,          tm115=0,    tm116=0.01399,
                tm121=0.00079260,    tm122=0.00001539,    tm123=0,          tm124=0,          tm125=0,    tm126=0.0009032,
                tm131=0,             tm132=0,             tm133=-0.08385,   tm134=-0.004435,  tm135=0,    tm136=0,
                tm141=0,             tm142=0,             tm143=-0.004435,  tm144=-0.0001854, tm145=0,    tm146=0,
                tm151=0,             tm152=0,             tm153=0,          tm154=0,          tm155=0,    tm156=0,
                tm161=0.01399,       tm162=0.0009032,     tm163=0,          tm164=0,          tm165=0,    tm166=-0.003659,

                tm211=0.09262,       tm212=0.009216,      tm213=0,          tm214=0,          tm215=0,    tm216=0.1261,
                tm221=0.009216,      tm222=-0.0001958,    tm223=0,          tm224=0,          tm225=0,    tm226=0.007306,
                tm231=0,             tm232=0,             tm233=-0.5359,    tm234=-0.01625,   tm235=0,    tm236=0,
                tm241=0,             tm242=0.0000,        tm243=-0.01625,   tm244=-0.0004059, tm245=0,    tm246=0,
                tm251=0,             tm252=0,             tm253=0,          tm254=0,          tm255=0,    tm256=0,
                tm261=0.1261,        tm262=0.007306,      tm263=0,          tm264=0,          tm265=0,    tm266=-0.03866,

                tm311=0,             tm312=0,             tm313=-0.04767,   tm314=-0.002516,  tm315=0,    tm316=0,
                tm321=0,             tm322=0,             tm323=-0.004032,  tm324=-0.0001378, tm325=0,    tm326=0,
                tm331=-0.04767,      tm332=-0.004032,     tm333=0,          tm334=0,          tm335=0,    tm336=-0.01546,
                tm341=-0.002516,     tm342=-0.0001378,    tm343=0,          tm344=0,          tm345=0,    tm346=-0.0006162,
                tm351=0,             tm352=0,             tm353=0,          tm354=0,          tm355=0,    tm356=0,
                tm361=0,             tm362=0,             tm363=-0.01546,   tm364=-0.0006162, tm365=0,    tm366=0,

                tm411=0,             tm412=0,             tm413=-0.8024,    tm414=-0.05174,   tm415=0,    tm416=0,
                tm421=0,             tm422=0,             tm423=-0.05203,   tm424=-0.002224,  tm425=0,    tm426=0,
                tm431=-0.8024,       tm432=-0.05203,      tm433=0,          tm434=0,          tm435=0,    tm436=-0.1824,
                tm441=-0.05174,      tm442=-0.002224,     tm443=0,          tm444=0,          tm445=0,    tm446=-0.007719,
                tm451=0,             tm452=0,             tm453=0,          tm454=0,          tm455=0,    tm456=0,
                tm461=0,             tm462=0,             tm463=-0.1824,    tm464=-0.007719,  tm465=0,    tm466=0
                '''

