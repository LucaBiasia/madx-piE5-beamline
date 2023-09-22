import numpy as np
from globals import *


# MY BEAMLINE
POS_IN = position['posQSF41'] - 0.250               # Initial position
POS_OUT = position['posMEGCOL'] - 0.371 + 0.195     # Final position
L_TOT = POS_OUT - POS_IN                            # Total length of the sequence
DX = 0.01                                           # slice length (m)

# DIRECTION
REVERSE = -1                                         # 1 from UPSTREAM to DOWNSTREAM (forward)    or    -1 from DOWNSTREAM to UPSTREAM (backwards)

# SAMPLING OF THE INITIAL TRACKS
N_SAMPLE = 1000                                     # Number of sampling (should be ~100 times greater than the desired number of tracks)
N_G4BL = 1                                          # Number of G4BL max number of particles (leave 1 when not using G4BL)

FIELDMAP_PATH = '/home/lucabiasia/Scrivania/piE5/fieldmaps'

# MY MOMENTUM
P0 = 0.028                                          # Nominal momentum [GeV]
DP0 = P0*0.0                                        # Nominal momentum spread [GeV]
f = 0.3/P0                                          # conversion factor from T/m to m^-2 (1/(B*rho))

# MY INITIAL PHASE SPACE
SIGMAX = 0.15755                                    # m
SIGMAXP = 0.073994                                  # rad
SIGMAY = 0.21099                                    # m
SIGMAYP = 0.120313                                  # rad
RHOX = 0.833                                        # x correlation
RHOY = 0.0001                                       # y correlation

ex = SIGMAX*SIGMAXP*np.sqrt(1-RHOX**2)              # x emittance
ey = SIGMAY*SIGMAYP*np.sqrt(1-RHOY**2)              # y emittance
alfax = -RHOX*SIGMAX*SIGMAXP/ex
alfay = -RHOY*SIGMAY*SIGMAYP/ey
betax = SIGMAX**2/ex
betay = SIGMAY**2/ey
gammax = (1 + alfax**2)/betax
gammay = (1 + alfay**2)/betay










################################################################## SUGGESTED INITIAL PHASE SPACES ###################################################################

# 25cm before QSF41 center---------------------------
# SIGMAX = 0.15755                                    # m
# SIGMAXP = 0.073994                                  # rad
# SIGMAY = 0.21099                                    # m
# SIGMAYP = 0.120313                                  # rad
# RHOX = 0.833                                        # x correlation
# RHOY = 0.0001 

# At Collimator---------------------------------------
# sigmax = 0.0268
# sigmaxp = 0.143983
# sigmay = 0.01115
# sigmayp = 0.037238
# RHOX = 0.018
# RHOY = 0.141

# AFTER AST--------------------------------------------
# sigmax = 0.015/2.35    # m
# sigmaxp = 0.450/2.35    # rad
# sigmay = 0.02/2.35    # m
# sigmayp = 0.120/2.35    # rad
# rhox = 0        # x correlation
# rhoy = 0        # y correlation

# AT PILL1---------------------------------------------
# sigmax = 0.0187
# sigmaxp = 0.0471
# sigmay = 0.0203
# sigmayp = 0.0349
# rhox = -0.014
# rhoy = 0.398

# 25cm before QSK41 center (2023)-----------------------
# sigmax = 0.0473
# sigmaxp = 0.0243
# sigmay = 0.0614
# sigmayp = 0.0174
# rhox = 0.567
# rhoy = 0.714 
# ex = 879*1e-6
# ey = 649*1e-6

# 25cm before QSK41 center (old)------------------------
# sigmax = 0.0448/2
# sigmaxp = 0.0205
# sigmay = 0.0504/2
# sigmayp = 0.0183
# rhox = 0.288
# rhoy = 0.710

# PHASE SPACE AT QSF41 CENTER----------------------------
# r12 = 0.833        # x correlation
# r34 = 0.0001      # y correlation
# chix = np.arcsin(r12)
# chiy = np.arcsin(r34)

# sigmax = 0.15755    # m
# sigmaxp = 0.073994*r12   # rad
# sigmay = 0.21099    # m
# sigmayp = 0.120313*r34   # rad

# ex = sigmax*sigmaxp*np.cos(chix)
# ey = sigmay*sigmayp*np.cos(chiy)
# alfax = -np.tan(chix)
# alfay = -np.tan(chiy)
# betax = sigmax**2/ex
# betay = sigmay**2/ey
# gammax = sigmaxp**2/ex
# gammay = sigmayp**2/ey

# rhox = np.tan(chix)*ex/(sigmax*sigmaxp)
# rhoy = -alfay*ey/(sigmay*sigmayp)




