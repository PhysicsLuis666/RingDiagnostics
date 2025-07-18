#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 14 14:05:59 2025

to produce x-x', y-y' and z-dE or x-y plots

@author: luisruiz
"""
import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def pathdirect(folder_path,filename,outputfolder):
    
    os.makedirs(outputfolder, exist_ok=True)
    
    output_path1 = os.path.join(outputfolder, f'Figure_{filename}.png')

    df = pd.read_csv(os.path.join(folder_path, filename), sep="\s+",comment='%', header=None, skip_blank_lines=True, usecols=range(6))
    
    X = df.to_numpy()
    X[:, 0:4] *= 1000  # scale x, x', y, y' to mm / mrad
    
    x, xp = X[:,0], X[:,1]
    y, yp = X[:,2], X[:,3]
    z, dE = X[:,4], X[:,5]
    
    x -= np.mean(x)
    xp -= np.mean(xp)
    y -= np.mean(y)
    yp -= np.mean(yp)
    z -= np.mean(z)
    dE -= np.mean(dE)
    
    #CovX = np.cov(X.T)
    #rmsvals = np.std(X, axis=0)
    
    #covxxp = CovX[0:2,0:2]
    ##covyyp = CovX[2:4,2:4]
    #covzdE = CovX[4:6,4:6]
    
    
    #------------------------------------------------------
    
    # In[9]:
    
    #setting up normalization
    #-----------------------------------------------------
    """
    e_x = np.sqrt(np.linalg.det(covxxp))
    e_y = np.sqrt(np.linalg.det(covyyp))
    e_z = np.sqrt(np.linalg.det(covzdE))
    
    rx = 2*rmsvals[0]
    rxp = 2*covxxp[0,1]/rmsvals[0]
    eps_x = 4*e_x
    
    ry = 2*rmsvals[2]
    ryp = 2*covyyp[0,1]/rmsvals[2]
    eps_y = 4*e_y
    
    rz = 2*rmsvals[4]
    rzp = 2*covzdE[0,1]/rmsvals[4]
    eps_z = 4*e_z
    
    
    #normalized coordinates 
    #--------------------------------------------------------
    # Now safe to normalize:
    xn = x / rx
    xpn = (rx*xp - rxp*x) / eps_x
    yn = y / ry
    ypn = (ry*yp - ryp*y) / eps_y

    
    #defining fraction of particles outside core
    #----------------------------------------------------------
    frac = .059  #float(input("fraction of particles outside of rms core?:"))
    rfrac = 1+frac
    rxy = xn**2 + yn**2
    rxxp = xn**2 + xpn**2
    ryyp = yn**2 + ypn**2
    rxpnypn = xpn**2 + ypn**2
    
    H_xy = np.sum((rxy > rfrac))/len(rxy)
    H_xxp = np.sum((rxxp > rfrac))/len(rxxp)
    H_yyp = np.sum((ryyp > rfrac))/len(ryyp)
    H_xpnypn = np.sum((rxpnypn > rfrac))/len(rxpnypn)
    
    print(rxpnypn)
    
    xn -= np.mean(xn)
    xpn -= np.mean(xpn)
    yn -= np.mean(yn)
    ypn -= np.mean(ypn)"""
    
    #----------------------------------------------------------------------
    #plot histograms plots
    #-----------------------------------------------------------------------
    """
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # x vs x'
    ax1.hist2d(x, xp, bins=128, cmap='Greys')
    ax1.set_xlim(-40, 40)
    ax1.set_ylim(-3.5, 3.5)
    ax1.set_xlabel("x mm")
    ax1.set_ylabel("x' mrad")
    ax1.set_title("x' vs x")
    #ax1.set_aspect('equal', adjustable='box')
    
    # x vs y
    ax2.hist2d(x, y, bins=128, cmap="Greys")
    ax2.set_xlim(-40, 40)
    ax2.set_ylim(-45, 45)
    ax2.set_xlabel("x mm")
    ax2.set_ylabel("y mm")
    ax2.set_title("y vs x")
    ax2.set_aspect('equal', adjustable='box')
    
    # y vs y'
    ax3.hist2d(y, yp, bins=128, cmap="Greys")
    ax3.set_xlim(-40, 40)
    ax3.set_ylim(-3.5, 3.5)
    ax3.set_xlabel("y mm")
    ax3.set_ylabel("y' mrad")
    ax3.set_title("y' vs y")
    #ax3.set_aspect('equal', adjustable='box')"""
    
    
    #----------------------------------------------------------------------
    #plot histograms plots
    #-----------------------------------------------------------------------
    """
    
    fig1 = plt.figure(figsize=(12,12))
    gs1 = gridspec.GridSpec(4,5)
     
    #x'
    ax1 = fig1.add_subplot(gs1[0,0])
    ax1.hist(x, bins = 128, color = 'Grey', histtype = 'step')
    ax1.set_xlim(-45, 45)
    #ax1.set_ylim(-4, 4)
    #ax1.tick_params(axis = 'x', pad = 20)
    #ax1.set_aspect('equal')
    ax1.set_xlabel("x mm")
    ax1.set_title("1d Hist x ")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.yaxis.tick_left()
    ax1.xaxis.tick_bottom()
    
    ax2 = fig1.add_subplot(gs1[1,0])
    ax2.hist2d(x, xp, bins = 128, cmap = 'Greys')
    ax2.set_xlim(-45, 45)
    ax2.set_ylim(-4, 4)
    #ax1.tick_params(axis = 'x', pad = 20)
    #ax1.set_aspect('equal')
    ax2.set_xlabel("x mm")
    ax2.set_ylabel("x' mrad")
    ax2.set_title("x' vs x")
    
    
    ax2 =fig1.add_subplot(gs1[1,1])
    ax2.hist(xp, bins=128, color = 'grey', histtype='step')
    ax2.set_xlim(-4,4)
    ax2.set_xlabel("x' mrad")
    ax2.set_title("1d Hist x'")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.yaxis.tick_left()
    ax2.xaxis.tick_bottom()
    
    #y
    ax3 = fig1.add_subplot(gs1[2,0])
    ax3.hist2d(x,y, bins = 128, cmap = 'Greys')
    ax3.set_xlim(-45, 45)
    ax3.set_ylim(-45, 45)
    ax3.set_xlabel("x mm")
    ax3.set_ylabel("y mm")
    ax3.set_title(" y vs x")
    
    ax3 = fig1.add_subplot(gs1[2,1])
    ax3.hist2d(xp,y, bins = 128, cmap = 'Greys')
    ax3.set_xlim(-4, 4)
    ax3.set_ylim(-45, 45)
    ax3.set_xlabel("x' mrad")
    ax3.set_ylabel("y mm")
    ax3.set_title(" y vs x'")
    
    ax3 = fig1.add_subplot(gs1[2,2])
    ax3.hist(y, bins = 128, color = 'grey', histtype='step')
    ax3.set_xlim(-45,45)
    ax3.set_xlabel('y mm')
    ax3.set_title("1d Hist y")
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.yaxis.tick_left()
    ax3.xaxis.tick_bottom()
    
        

    
 
    
    
    # y'
    
    ax4 = fig1.add_subplot(gs1[3,0])
    ax4.hist2d(x,yp, bins = 128, cmap = 'Greys')
    ax4.set_xlim(-45, 45)
    ax4.set_ylim(-4, 4)
    ax4.set_xlabel("x mm")
    ax4.set_ylabel("y' mrad")
    ax4.set_title("y' vs x")
    
    ax4 = fig1.add_subplot(gs1[3,1])
    ax4.hist2d(xp,yp, bins = 128, cmap = 'Greys')
    ax4.set_xlim(-4, 4)
    ax4.set_ylim(-4, 4)
    ax4.set_xlabel("x' mrad")
    ax4.set_ylabel("y' mrad")
    ax4.set_title("y' vs x' ")
    
    ax4 = fig1.add_subplot(gs1[3,2])
    ax4.hist2d(y,yp, bins = 128, cmap = 'Greys')
    ax4.set_xlim(-45, 45)
    ax4.set_ylim(-4, 4)
    ax4.set_xlabel("y mm")
    ax4.set_ylabel("y' mrad")
    ax4.set_title("y' v y ")
    
    ax4 = fig1.add_subplot(gs1[3,3])
    ax4.hist(y, bins = 128, color = 'Grey', histtype='step')
    ax4.set_xlim(-45,45)
    ax4.set_xlabel("y' mrad")
    ax4.set_title("1d Hist y'")
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    ax4.yaxis.tick_left()
    ax4.xaxis.tick_bottom()
    
    #ax5 = fig1.add_subplot(gs1[3,4])
    #ax5.hist2d(z,dE, bins = 128, cmap = 'Greys')
    #ax4.set_xlim(-4, 4)
    #ax4.set_ylim(-45, 45)
    #ax5.set_xlabel("z m")
    #ax5.set_ylabel("dE GeV")
    #ax5.set_title("dE vs z'")

    fig1.suptitle(f'Phase space of {filename}', fontsize=12)
    fig1.tight_layout()
    fig1.savefig(output_path1)
    plt.close(fig1)
    """
    fig1 = plt.figure(figsize=(10, 10))  # Wider canvas
    gs1 = gridspec.GridSpec(4, 4)
    
    # x 1D
    ax_x = fig1.add_subplot(gs1[0, 0])
    ax_x.hist(x, bins=128, color='grey', histtype='step')
    ax_x.set_xlim(-50, 50)
    ax_x.set_xlabel("x mm")
    ax_x.set_title("1D Hist x")
    ax_x.spines["top"].set_visible(False)
    ax_x.spines["right"].set_visible(False)
  
    
    # x-x'
    ax_xxp = fig1.add_subplot(gs1[1, 0])
    ax_xxp.hist2d(x, xp, bins=128, cmap='Greys')
    ax_xxp.set_xlim(-50,50)
    ax_xxp.set_ylim(-5, 5)
    ax_xxp.set_xticks(np.linspace(-50,50, 5))  # for example: -45, -30, ..., 45
    ax_xxp.set_yticks(np.linspace(-5, 5, 5))   
    ax_xxp.set_xlabel("x mm")
    ax_xxp.set_ylabel("x' mrad")
    ax_xxp.set_title("x' vs x")
    
    
    # x' 1D
    ax_xp = fig1.add_subplot(gs1[1, 1])
    ax_xp.hist(xp, bins=128, color='grey', histtype='step')
    ax_xp.set_xlim(-5, 5)
    ax_xp.set_xlabel("x' mrad")
    ax_xp.set_title("1D Hist x'")
    ax_xp.spines["top"].set_visible(False)
    ax_xp.spines["right"].set_visible(False)

    
    # x-y
    ax_xy = fig1.add_subplot(gs1[2, 0])
    ax_xy.hist2d(x, y, bins=128, cmap='Greys')
    ax_xy.set_xlim(-50, 50)
    ax_xy.set_ylim(-50, 50)
    ax_xy.set_xticks(np.linspace(-50, 50, 5))  # for example: -45, -30, ..., 45
    ax_xy.set_yticks(np.linspace(-50, 50, 5))  
    ax_xy.set_xlabel("x mm")
    ax_xy.set_ylabel("y mm")
    ax_xy.set_title("y vs x")
    
    # x'-y
    ax_xpy = fig1.add_subplot(gs1[2, 1])
    ax_xpy.hist2d(xp, y, bins=128, cmap='Greys')
    ax_xpy.set_xlim(-5, 5)
    ax_xpy.set_ylim(-50,50)
    ax_xpy.set_yticks(np.linspace(-50, 50, 5)) 
    ax_xpy.set_xticks(np.linspace(-5, 5, 5))  # for example: -45, -30, ..., 45
    #ax_xpy.set_yticks(np.linspace(-45, 45, 5))  
    ax_xpy.set_xlabel("x' mrad")
    ax_xpy.set_ylabel("y mm")
    ax_xpy.set_title("y vs x'")
    
    
    # y 1D
    ax_y = fig1.add_subplot(gs1[2, 2])
    ax_y.hist(y, bins=128, color='grey', histtype='step')
    ax_y.set_xlim(-50, 50)
    ax_y.set_xlabel("y mm")
    ax_y.set_title("1D Hist y")
    ax_y.spines["top"].set_visible(False)
    ax_y.spines["right"].set_visible(False)
    
    # y'-x
    ax_ypx = fig1.add_subplot(gs1[3, 0])
    ax_ypx.hist2d(x, yp, bins=128, cmap='Greys')
    ax_ypx.set_xlim(-50, 50)
    ax_ypx.set_ylim(-5, 5)
    ax_ypx.set_xticks(np.linspace(-50, 50, 5)) 
    ax_ypx.set_yticks(np.linspace(-5, 5, 5)) 
    ax_ypx.set_xlabel("x mm")
    ax_ypx.set_ylabel("y' mrad")
    ax_ypx.set_title("y' vs x")
    
    # y'-x'
    ax_ypxp = fig1.add_subplot(gs1[3, 1])
    ax_ypxp.hist2d(xp, yp, bins=128, cmap='Greys')
    ax_ypxp.set_xlim(-5, 5)
    ax_ypxp.set_ylim(-5, 5)
    ax_ypxp.set_xticks(np.linspace(-5, 5, 5)) 
    ax_ypxp.set_yticks(np.linspace(-5, 5, 5)) 
    ax_ypxp.set_xlabel("x' mrad")
    ax_ypxp.set_ylabel("y' mrad")
    ax_ypxp.set_title("y' vs x'")
    
    # y'-y
    ax_ypy = fig1.add_subplot(gs1[3, 2])
    ax_ypy.hist2d(y, yp, bins=128, cmap='Greys')
    ax_ypy.set_xlim(-50, 50)
    ax_ypy.set_ylim(-5, 5)
    ax_ypy.set_xticks(np.linspace(-50, 50, 5)) 
    ax_ypy.set_yticks(np.linspace(-5, 5, 5)) 
    ax_ypy.set_xlabel("y mm")
    ax_ypy.set_ylabel("y' mrad")
    ax_ypy.set_title("y' vs y")
    
    # y' 1D
    ax_yp = fig1.add_subplot(gs1[3, 3])
    ax_yp.hist(yp, bins=128, color='grey', histtype='step')
    ax_yp.set_xlim(-5, 5)
    ax_yp.set_xlabel("y' mrad")
    ax_yp.set_title("1D Hist y'")
    ax_yp.spines["top"].set_visible(False)
    ax_yp.spines["right"].set_visible(False)
    
    # Title & Save
    fig1.suptitle('Phase Space', fontsize=14)
    fig1.tight_layout(pad=2.0)
    fig1.savefig(output_path1)
    plt.close(fig1)


    
    

turns = 980 #int(input("how many turns?"))
frequency = 250 #int(input("what was the write frequency?"))

for i in list(range(0,turns, frequency)) + [turns-1]:
    filename = f"bunch_{i:04d}"
    print(f"Processing {filename}")
    pathdirect('outputs2-2/nullwaveform1',filename,
               '/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs2-2/nullwaveform1/Analysis1-1')


