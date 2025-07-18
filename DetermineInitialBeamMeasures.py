#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 14:40:18 2025

@author: luisruiz
"""

import math
import os
import sys
import time
import json

import numpy as np
import matplotlib.pyplot as plt

from orbit.core import orbit_mpi
from orbit.core.bunch import Bunch
from orbit.core.bunch import BunchTwissAnalysis
from orbit.core.spacecharge import Boundary2D
from orbit.core.spacecharge import LSpaceChargeCalc
from orbit.core.spacecharge import SpaceChargeCalc2p5D
from orbit.core.spacecharge import SpaceChargeCalcSliceBySlice2D
from orbit.aperture import TeapotApertureNode
from orbit.aperture import CircleApertureNode
from orbit.aperture import RectangleApertureNode
from orbit.aperture import EllipseApertureNode
from orbit.bumps import TeapotBumpNode
from orbit.bumps import TDTeapotSimpleBumpNode
from orbit.collimation import TeapotCollimatorNode
from orbit.collimation import addTeapotCollimatorNode
from orbit.foils import TeapotFoilNode
from orbit.impedances import addImpedanceNode
from orbit.impedances import LImpedance_Node
from orbit.impedances import TImpedance_Node
from orbit.injection import TeapotInjectionNode
from orbit.injection import JohoTransverse
from orbit.injection import JohoLongitudinal
from orbit.injection import SNSESpreadDist
from orbit.kickernodes import flatTopWaveform
from orbit.kickernodes import SquareRootWaveform
from orbit.lattice import AccLattice
from orbit.lattice import AccNode
from orbit.lattice import AccActionsContainer
from orbit.rf_cavities import RFNode
from orbit.rf_cavities import addRFNode
from orbit.space_charge.sc1d import SC1D_AccNode
from orbit.space_charge.sc1d import addLongitudinalSpaceChargeNode
from orbit.space_charge.sc2p5d import SC2p5D_AccNode
from orbit.space_charge.sc2p5d import setSC2p5DAccNodes
from orbit.space_charge.sc2dslicebyslice import SC2DSliceBySlice_AccNode
from orbit.space_charge.sc2dslicebyslice import setSC2DSliceBySliceAccNodes
from orbit.teapot import TEAPOT_Ring
from orbit.teapot import DriftTEAPOT
from orbit.utils.consts import mass_proton
from orbit.utils.consts import speed_of_light

# local



#begin code here 
#cov matrix
#----------------------------------------------------------------------------------------
def get_bunch_cov(                              
    bunch: Bunch, dispersion_flag: bool = False, emit_norm_flag: bool = False
) -> np.ndarray:                              #this is a function that obtains the covaraint matrix found in the github ripository
    order = 2
    twiss_calc = BunchTwissAnalysis()
    twiss_calc.computeBunchMoments(
        bunch,
        order,
        int(dispersion_flag),
        int(emit_norm_flag),
    )

    S = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1):
            S[i, j] = twiss_calc.getCorrelation(j, i)
            S[j, i] = S[i, j]
    return S
#------------------------------------------------------------------------------

def get_bunch_coords(bunch: Bunch):          #defines the function that obtainst the particles coordinates
    X = np.zeros((bunch.getSize(),6))
    for i in range(bunch.getSize()):
        x = bunch.x(i)
        y = bunch.y(i)
        z = bunch.z(i)
        xp = bunch.xp(i)
        yp = bunch.yp(i)
        de = bunch.dE(i)
        X[i, :] = (x,xp,y,yp,z,de)
    return X   

#------------------------------------------------------------------------------

power = 1.7  # [MW]
energy = 1.3  # [GeV]
frequency = 60.0  # [Hz]

intensity = 1.38E14 # previous was a bit different refere to track.py, 1.38E14, 2.0E14 , 1.38E14, 1.38e15

minipulses_per_pulse = 980  #1044       
minipulse_intensity = intensity / minipulses_per_pulse
macrosize = minipulse_intensity / 500

bunch = Bunch()
bunch.mass(mass_proton)
bunch.macroSize(macrosize)
bunch.getSyncParticle().kinEnergy(1.3)  # [GeV]


# Initialize lattice
# --------------------------------------------------------------------------------------

lattice = TEAPOT_Ring()
lattice.readMAD("inputs/sns_ring_mad.lattice", "RINGINJ")
lattice.initialize()



nodes = lattice.getNodes()

for node in nodes:
    print(node.getName(), node)


# Set injection kicker waveforms
# --------------------------------------------------------------------------------------

# Collect kicker nodes
hkick10_node = lattice.getNodeForName("IKICKH_A10")
vkick10_node = lattice.getNodeForName("IKICKV_A10")
hkick11_node = lattice.getNodeForName("IKICKH_A11")
vkick11_node = lattice.getNodeForName("IKICKV_A11")
hkick12_node = lattice.getNodeForName("IKICKH_A12")
vkick12_node = lattice.getNodeForName("IKICKV_A12")
hkick13_node = lattice.getNodeForName("IKICKH_A13")
vkick13_node = lattice.getNodeForName("IKICKV_A13")


#Limit Strengths
limit_hkicker10 = limit_hkicker13 = 12.26e-03
limit_vkicker10 = limit_vkicker13 = 12.25e-03
limit_hkicker11 = limit_hkicker12 = 6.80e-03
limit_vkicker11 = limit_vkicker12 = 6.79e-03

bump_hkicker10 =0.9956521739130436  
bump_vkicker10 =0.6278260869565218
bump_hkicker11 =0.5469565217391305
bump_vkicker11 =0.6474890272627937 
bump_hkicker12 =0.6478260869565218    
bump_vkicker12 =0.6817391304347827 
bump_hkicker13 =.99 #1
bump_vkicker13 =0.6396521739130435 

#Horizontal ~ si=0.99 sf=0.48
#Vertical si=0.58 sf=0.29
# Set kicker strengths
strength_hkicker10 = limit_hkicker10 * bump_hkicker10
strength_vkicker10 = limit_vkicker10 * bump_vkicker10
strength_hkicker11 = -limit_hkicker11 * bump_hkicker11
strength_vkicker11 = -limit_vkicker11 * bump_vkicker11
strength_hkicker12 = -limit_hkicker12 * bump_hkicker12
strength_vkicker12 = -limit_vkicker12 * bump_vkicker12
strength_hkicker13 = limit_hkicker13 * bump_hkicker13
strength_vkicker13 = limit_vkicker13 * bump_vkicker13

hkick10_node.setParam("kx", strength_hkicker10)
vkick10_node.setParam("ky", strength_vkicker10)
hkick11_node.setParam("kx", strength_hkicker11)
vkick11_node.setParam("ky", strength_vkicker11)
hkick12_node.setParam("kx", strength_hkicker12)
vkick12_node.setParam("ky", strength_vkicker12)
hkick13_node.setParam("kx", strength_hkicker13)
vkick13_node.setParam("ky", strength_vkicker13)

# Set kicker waveforms
tih = -0.001
tiv = -0.002
tf =  0.001
si =  1
sih = .99
siv = .77  #.77, .58
sfh =  0.5538  #378,.457 .48, 0.557573720104 , 0.556, .557, .554 
sfv =  0.423 #.29,.359,.406 #.38, .42, .304567 .432659196502 , 0.432659196502 , .43, 0.425212690535, .422

sync_part = bunch.getSyncParticle()
hkickerwave = SquareRootWaveform(sync_part, lattice.getLength(), tih, tf, sih, sfh)
vkickerwave = SquareRootWaveform(sync_part, lattice.getLength(), tiv, tf, siv, sfv)

hkick10_node.setWaveform(hkickerwave)
vkick10_node.setWaveform(vkickerwave)
hkick11_node.setWaveform(hkickerwave)
vkick11_node.setWaveform(vkickerwave)
hkick12_node.setWaveform(hkickerwave)
vkick12_node.setWaveform(vkickerwave)
hkick13_node.setWaveform(hkickerwave)
vkick13_node.setWaveform(vkickerwave)



# Define minipulse distribution function
# --------------------------------------------------------------------------------------

lattice_length = lattice.getLength()
sync_part = bunch.getSyncParticle()

order = 9.0
alphax = 0.064
betax = 10.056
alphay = 0.063
betay = 10.815
emitlim = 0.221 * 2 * (order + 1) * 1.0e-6
xcenterpos = 0.0486
xcentermom = 0.0
ycenterpos = 0.046
ycentermom = 0.0

zlim = 139.68 * lattice_length / 360.0
zmin = -zlim
zmax = zlim
tailfraction = 0.0
emean = sync_part.kinEnergy()
esigma = 0.0005
etrunc = 1.0
emin = sync_part.kinEnergy() - 0.0025
emax = sync_part.kinEnergy() + 0.0025
ecmean = 0.0
ecsigma = 0.000000001
ectrunc = 1.0
ecmin = -0.0035
ecmax = 0.0035
ecdrifti = 0.0
ecdriftf = 0.0
time_per_turn = lattice_length / (sync_part.beta() * speed_of_light)
drifttime = 1000.0 * 981 * time_per_turn
ecparams = (ecmean, ecsigma, ectrunc, ecmin, ecmax, ecdrifti, ecdriftf, drifttime)
esnu = 100.0
esphase = 0.0
esmax = 0.0
nulltime = 0.0
esparams = (esnu, esphase, esmax, nulltime)

inj_dist_x = JohoTransverse(order, alphax, betax, emitlim, xcenterpos, xcentermom)
inj_dist_y = JohoTransverse(order, alphay, betay, emitlim, ycenterpos, ycentermom)
inj_dist_z = SNSESpreadDist(lattice_length, zmin, zmax, tailfraction, sync_part, emean, esigma, etrunc, emin, emax, ecparams, esparams)
X, Xp, Y, Yp, Z, dE = [], [], [], [], [], []

for i in range(100_0000):
    xpos, xmom = inj_dist_x.getCoordinates()
    ypos, ymom = inj_dist_y.getCoordinates()
    zinj, dEinj = inj_dist_z.getCoordinates()

    X.append(xpos)
    Xp.append(xmom)
    Y.append(ypos)
    Yp.append(ymom)
    Z.append(zinj)
    dE.append(dEinj)
    

# Convert to numpy arrays
X = np.array(X)
Xp = np.array(Xp)
Y = np.array(Y)
Yp = np.array(Yp)
Z = np.array(Z)
dE = np.array(dE)

# Stack as columns: shape will be (100000, 6)
phase_space_array = np.column_stack((X, Xp, Y, Yp, Z, dE))

# Access columns (not rows)
x  = phase_space_array[:, 0]
xp = phase_space_array[:, 1]
y  = phase_space_array[:, 2]
yp = phase_space_array[:, 3]
z  = phase_space_array[:, 4]
de = phase_space_array[:, 5]

# Covariance matrix (6x6)
CovM = np.cov(phase_space_array.T)

#print(CovM)
print(f"this is the exrms from the pyorbit method: ex = {inj_dist_x.emitrms}")
print(f"this is the eyrms from the pyorbit method: ey = {inj_dist_y.emitrms}")

covX = CovM[0:2,0:2]
covY = CovM[2:4,2:4]
covZ = CovM[4:6,4:6]

print(f"this is matrix x:{covX}")
print(f"this is matrix y:{covY}")
print(f"this is matriz z:{covZ}")

eps_xrms = np.sqrt(np.linalg.det(covX))
eps_yrms = np.sqrt(np.linalg.det(covY))
eps_zrms = np.sqrt(np.linalg.det(covZ))

print(f"rms emittances ex = {eps_xrms*1000000}, ey = {eps_yrms*1000000}, ez = {eps_zrms}")
rmsvals = np.std(phase_space_array, axis=0)
#xmean = np.sum(x)/len(x)
#xrms = np.sqrt(np.sum((x-xmean)**2)/(len(x)-1))
varx = covX[0,0]  #<x^2>
varxxp = covX[0,1] #<xx'>
varxp = covX[1,1] #<x'^2>
betax = varx/eps_xrms
vary = covY[0,0]  #<x^2>
varyyp = covY[0,1] #<xx'>
varyp = covY[1,1] #<x'^2>
betay = vary/eps_yrms
xprms = np.sqrt(varxp)
xrms = np.sqrt(varx)
yrms = np.sqrt(vary)
#emit_xrms = np.sqrt(varx*varxp - varxxp**2)
#ymean = np.sum(y)/len(y)
#yrms = np.sqrt(np.sum((y-ymean)**2)/(len(y)-1))
#yrms = np.sqrt(covY[0,0])
#print(f"calculating directly xrms:{xrms}")
#print(f"using cov matrix:{xrms2}")
#print(f"rms y direct:{yrms}")
#print(f"y rms using cov matrix {yrms2}")
vardE = covZ[1,1]  #< dE ^2>
varz = covZ[0,0]  # <z^2>
varzdE =covZ[0,1]
rmsdE = np.sqrt(vardE)
rmsz = np.sqrt(varz)
sigmax = np.sqrt(betax*eps_xrms)
sigmay = np.sqrt(betay*eps_yrms)
rmsenergyspread = rmsdE/1.3
rx = 2*xrms*1000
ry = 2*yrms*1000

print(f"rms dE:{rmsdE}")
print(f"beta x = {betax}")
print(f"beta y = {betay}")
print(f"beamsizex = {sigmax*1000}, beamsizey = {sigmay*1000}")
print(f"rms energy spread:{rmsenergyspread}")
print(f"rms z = {rmsz}")
print(f"Rms edge radius x = {rx}")
print(f"Rms edge radius y = {ry}")

"""
for i in range(10):
    
    lattice.trackBunch(bunch)
    
    X = get_bunch_coords(bunch)
    X = X * 1000.0

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 10))
    
    ax1.scatter(X[:, 0], X[:, 1], s=1, color='blue')
    ax1.set_xlabel("x [mm]")
    ax1.set_ylabel("x' [mrad]")
    
    ax2.scatter(X[:, 2], X[:, 3], s=1, color='red')
    ax2.set_xlabel("y [mm]")
    ax2.set_ylabel("y' [mrad]")
    
    ax3.scatter(X[:, 0], X[:, 2], s=1, color='red')
    ax3.set_xlabel("x [mm]")
    ax3.set_ylabel("y [mrad]")
    
    plt.tight_layout()
    plt.show()
    
    plt.tight_layout()
    plt.show()
    """



