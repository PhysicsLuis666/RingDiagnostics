#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 10:08:39 2025

@author: luisruiz
"""

"""
New Analysis Script with multiple sets of data

"""


import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def process_data(file_path):
    
    df = pd.read_csv(file_path, sep="\s+",comment='%', header=None, skip_blank_lines=True, usecols=range(6))
    
    X = df.to_numpy()
    
    X[:, 0:4] *= 1000  # scale x, x', y, y' to mm / mrad
    
    X -= np.mean(X, axis=0)
    
    return X
    
    
    

def comparison_plot(X1,X2,X3,X4,X5,X6,X7, filename, output_dir): 
    x1, xp1, y1, yp1 = X1[:, 0], X1[:, 1], X1[:, 2], X1[:, 3] #FR1
    x2, xp2, y2, yp2 = X2[:, 0], X2[:, 1], X2[:, 2], X2[:, 3] #FR2
    x3, xp3, y3, yp3 = X3[:, 0], X3[:, 1], X3[:, 2], X3[:, 3] #FR3
    x4, xp4, y4, yp4 = X4[:, 0], X4[:, 1], X4[:, 2], X4[:, 3] #FR4
    x5, xp5, y5, yp5 = X5[:, 0], X5[:, 1], X5[:, 2], X5[:, 3] #FR5
    x6, xp6, y6, yp6 = X6[:, 0], X6[:, 1], X6[:, 2], X6[:, 3] #FR6
    x7, xp7, y7, yp7 = X7[:, 0], X7[:, 1], X7[:, 2], X7[:, 3] #FR7
    
    
    #fig, axs = plt.subplots(1, 3, figsize=(12, 12))
    """
    axs[0].hist(x1, bins=128, color='blue', density=False,alpha=0.5, label='SFH=.56,SFV=.418 simulation', histtype='step')   #log=True
    axs[0].hist(x2, bins=128, color='red', density=False, alpha=0.5, label='SFH=.56, SFV=.42 simulation', histtype='step')
    axs[0].set_title('x distribution comparison between two settings')
    axs[0].set_xlabel('x [mm]')
    axs[0].legend(fontsize=8,            # Smaller font
                  loc='upper right')
    
    axs[1].hist(x3, bins=128, color='blue', density=False,alpha=0.5, label='SFH=.554,SFV=.42 simulation', histtype='step')
    axs[1].hist(x4, bins=128, color='red', density=False, alpha=0.5, label='SFH=.555, SFV=.415 simulation', histtype='step')
    axs[1].set_title('x distribution comparison between two settings')
    axs[1].set_xlabel('x [mm]')
    axs[1].legend(fontsize=8,            # Smaller font
                  loc='upper right')
    
    axs[2].hist(x1, bins=128, color='blue', density=False,alpha=0.5, label='SFH=.56,SFV=.418 emittance values obtained', histtype='step')
    axs[2].hist(x4, bins=128, color='red', density=False, alpha=0.5, label='SFH=.555, SFV=.415 simulation emittance values obtained', histtype='step')
    axs[2].set_title('x distribution comparison between two settings')
    axs[2].set_ylabel("counts")
    axs[2].set_xlabel('x [mm]')
    axs[2].legend(fontsize=8,            # Smaller font
                  loc='upper right')     # Legend position
                  #bbox_to_anchor=(1, 1), # Offset position (if needed)
                  #ncol=1,                # One column
                  #frameon=True  )
    
      
    plt.suptitle(f'Comparison: {filename}')
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'compare_{filename}.png'))
    plt.close()"""
    
    
    fig1, axs1 = plt.subplots(2,2, figsize = (13,13))
    axs1[0,0].hist(x5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17 control room values obtained in simulation', histtype='step', log=False)
    axs1[0,0].hist(x1, bins=128, color='red', density=True, alpha=0.5, label='SFH=.56, SFV=.418 simulation', histtype='step', log=False)
    axs1[0,0].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    axs1[0,0].set_xlabel('x [mm]')
    axs1[0,0].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs1[1,0].hist(y5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17 control room values obtained in simulation', histtype='step', log=True)
    axs1[1,0].hist(y1, bins=128, color='red', density=True, alpha=0.5, label='SFH=.56, SFV=.418 simulation', histtype='step', log=True)
    axs1[1,0].set_title('y distribution comparison between two settings')
    #axs.set_xlabel()
    axs1[1,0].set_xlabel('y [mm]')
    axs1[1,0].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs1[0,1].hist(x5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17 control room values obtained in simulation', histtype='step', log=True)
    axs1[0,1].hist(x4, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.415 simulation', histtype='step', log=True)
    axs1[0,1].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    axs1[0,1].set_xlabel('x [mm]')
    axs1[0,1].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs1[1,1].hist(y5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17 control room values obtained in simulation', histtype='step', log=True)
    axs1[1,1].hist(y4, bins=128, color='red', density=True, alpha=0.5, label='SFH=.555, SFV=.415 simulation', histtype='step', log=True)
    axs1[1,1].set_title('xydistribution comparison between two settings')
    #axs.set_xlabel()
    axs1[1,1].set_xlabel('y [mm]')
    axs1[1,1].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    plt.suptitle(f'Comparison of control kicker strength values to simulaton values: {filename}')
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'compare control room values to simulation values_{filename}.png'))
    plt.close(fig1)
    
    """ Final Run 10 and coompare to run 5
    # "sih": 0.99,
    #"siv": 0.99,
    #"sfh": 0.56,
    #"sfv": 0.39,"""
    fig2, axs2 = plt.subplots(2,2, figsize = (13,13))
    axs2[0,0].hist(x5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17, SIV=.77 control room values obtained in simulation', histtype='step', log=False)
    axs2[0,0].hist(x7, bins=128, color='red', density=True, alpha=0.5, label='SFH=.56, SFV=.39, SIV=.99 simulation', histtype='step', log=False)
    axs2[0,0].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    axs2[0,0].set_xlabel('x [mm]')
    axs2[0,0].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs2[1,0].hist(y5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17m SIV=.77 control room values obtained in simulation', histtype='step', log=False)
    axs2[1,0].hist(y7, bins=128, color='red', density=True, alpha=0.5, label='SFH=.56, SFV=.39, =.99 ', histtype='step', log=False)
    axs2[1,0].set_title('y distribution comparison between two settings')
    #axs.set_xlabel()
    axs2[1,0].set_xlabel('y [mm]')
    axs2[1,0].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs2[0,1].hist(x1, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.56,SFV=.418, SIV=.77', histtype='step', log=False)
    axs2[0,1].hist(x7, bins=128, color='red', density=True, alpha=0.5, label='SFH=.56, SFV=.39 , SIV=.99', histtype='step', log=False)
    axs2[0,1].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    axs2[0,1].set_xlabel('x [mm]')
    axs2[0,1].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs2[1,1].hist(y1, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.56,SFV=.418, SIV=.77 ', histtype='step', log=False)
    axs2[1,1].hist(y7, bins=128, color='red', density=True, alpha=0.5, label='SFH=.56, SFV=.39, SIV=.99', histtype='step', log=False)
    axs2[1,1].set_title('ydistribution comparison between two settings')
    #axs.set_xlabel()
    axs2[1,1].set_xlabel('y [mm]')
    axs2[1,1].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    plt.suptitle(f'Comparison Vlaues: {filename}')
    plt.tight_layout()
    os.makedirs(output_dir1, exist_ok=True)
    plt.savefig(os.path.join(output_dir1, f'ComparingLargeSIVpoint.99-1_{filename}.png'))
    plt.close(fig2)
    
    """
    sih": 0.99,
    "siv": 1.25,
    "sfh": 0.55,
    "sfv": 0.36,"""
    
    fig3, axs3 = plt.subplots(2,2, figsize = (13,13))
    axs3[0,0].hist(x5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17, SIV=.77 control room values obtained in simulation', histtype='step', log=False)
    axs3[0,0].hist(x6, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.36, SIV=1.25 simulation', histtype='step', log=False)
    axs3[0,0].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    axs3[0,0].set_xlabel('x [mm]')
    axs3[0,0].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs3[1,0].hist(y5, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.48,SFV=.17m SIV=.77 control room values obtained in simulation', histtype='step', log=False)
    axs3[1,0].hist(y6, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.36, SIV=1.25 simulation', histtype='step', log=False)
    axs3[1,0].set_title('y distribution comparison between two settings')
    #axs.set_xlabel()
    axs3[1,0].set_xlabel('y [mm]')
    axs3[1,0].legend(fontsize=8,            # Smaller font
                  loc='lower center')
    
    axs3[0,1].hist(x1, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.56,SFV=.418, SIV=.77simulation', histtype='step', log=False)
    axs3[0,1].hist(x6, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.36 , SIV=1.25, simulation', histtype='step', log=False)
    axs3[0,1].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    axs3[0,1].set_xlabel('x [mm]')
    axs3[0,1].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    axs3[1,1].hist(y1, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.56,SFV=.418, SIV=.77 control room values obtained in simulation', histtype='step', log=False)
    axs3[1,1].hist(y6, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.36 simulation, SIV=1.25', histtype='step', log=False)
    axs3[1,1].set_title('y distribution comparison between two settings')
    #axs.set_xlabel()
    axs3[1,1].set_xlabel('y [mm]')
    axs3[1,1].legend(fontsize=6,            # Smaller font
                  loc='lower center')
    
    plt.suptitle(f'Comparison Vlaues2: {filename}')
    plt.tight_layout()
    os.makedirs(output_dir2, exist_ok=True)
    plt.savefig(os.path.join(output_dir2, f'ComparingLargeSIVpoint1.25_{filename}.png'))
    plt.close(fig3)
    
        
    fig4, axs4 = plt.subplots(1, 2, figsize=(10, 5))
    
    axs4[0].hist(x1, bins=128, color='blue', density=True, alpha=0.5,
                 label='SFH=.56, SFV=.418, SIV=.77', histtype='step', log=False)
    axs4[0].hist(x7, bins=128, color='red', density=True, alpha=0.5,
                 label='SFH=.56, SFV=.39, SIV=.99', histtype='step', log=False)
    axs4[0].set_title('x 1D Histogram')
    axs4[0].set_xlabel('x [mm]')
    axs4[0].set_ylabel("Probability Density")
    axs4[0].legend(fontsize=6, loc='lower center')
    
    axs4[1].hist(y1, bins=128, color='blue', density=True, alpha=0.5,
                 label='SFH=.56, SFV=.418, SIV=.77', histtype='step', log=False)
    axs4[1].hist(y7, bins=128, color='red', density=True, alpha=0.5,
                 label='SFH=.56, SFV=.39, SIV=.99', histtype='step', log=False)
    axs4[1].set_title('y 1D Histogram')
    axs4[1].set_xlabel('y [mm]')
    axs4[1].set_ylabel("Probability Density")
    axs4[1].legend(fontsize=6, loc='lower center')
    """
    axs4[2].hist(x1, bins=128, color='blue', density=False, alpha=0.5,
                 label='SFH=.56, SFV=.418, SIV=.77', histtype='step', log=False)
    axs4[2].hist(x7, bins=128, color='red', density=False, alpha=0.5,
                 label='SFH=.56, SFV=.39, SIV=.99', histtype='step', log=False)
    axs4[2].set_title('x distribution')
    axs4[2].set_xlabel('x [mm]')
    axs4[2].legend(fontsize=6, loc='lower center')
    
    axs4[3].hist(y1, bins=128, color='blue', density=False, alpha=0.5,
                 label='SFH=.56, SFV=.418, SIV=.77', histtype='step', log=False)
    axs4[3].hist(y7, bins=128, color='red', density=False, alpha=0.5,
                 label='SFH=.56, SFV=.39, SIV=.99', histtype='step', log=False)
    axs4[3].set_title('y distribution')
    axs4[3].set_xlabel('y [mm]')
    axs4[3].legend(fontsize=6, loc='lower center')"""
    
    plt.tight_layout()
        
    #axs4[0,1].hist(x4, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.56,SFV=.418, SIV=.77simulation', histtype='step', log=True)
    #axs4[0,1].hist(x6, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.36 , SIV=1.25, simulation', histtype='step', log=True)
   # axs4[0,1].set_title('x distribution comparison between two settings')
    #axs.set_xlabel()
    #axs4[0,1].set_xlabel('x [mm]')
    #axs4[0,1].legend(fontsize=6,            # Smaller font
                  #loc='lower center')
    
    #axs4[1,1].hist(y4, bins=128, color='blue', density=True,alpha=0.5, label='SFH=.555, SFV=.415 , SIV=.77', histtype='step', log=True)
    #axs4[1,1].hist(y6, bins=128, color='red', density=True, alpha=0.5, label='SFH=.55, SFV=.36, SIV=1.25', histtype='step', log=True)
    #axs4[1,1].set_title('xydistribution comparison between two settings')
    #axs.set_xlabel()
    #axs4[1,1].set_xlabel('y [mm]')
   # axs4[1,1].legend(fontsize=6,            # Smaller font
                  #loc='lower center')
    
    plt.suptitle(f"1D Histogram Comparisons at Different Kicker Settings")
    plt.tight_layout()
    os.makedirs(output_dir3, exist_ok=True)
    plt.savefig(os.path.join(output_dir3, f'ComparingFilesLargeSIVpoint.99toLowerSIVpoints_{filename}.png'))
    plt.close(fig4)
    
    
    
        
        
    
      

folder1 = "outputs2/FinalRun1"  #older run .418 sfv and .77 siv
folder2 = "outputs2/FinalRun2"
folder3 = "outputs2/FinalRun3"
folder4 = "outputs2/FinalRun4"
folder5 = "outputs2/FinalRun5"  #this is the one that Nick use after tuning with .77 siv and .17 sfv
folder11 = "outputs2/FinalRun11"
folder10 = "outputs2/FinalRun10" #FinalRun10 
output_dir1 = "/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs2/Analysis1-3"
output_dir2 = "/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs2/Analysis1-4"
output_dir3 = "/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs2/Analysis1-5"


files1 = set(os.listdir(folder1))
files2 = set(os.listdir(folder2))
files3 = set(os.listdir(folder3))
files4 = set(os.listdir(folder4))
files5 = set(os.listdir(folder5))
files6 = set(os.listdir(folder11))
files7 = set(os.listdir(folder10)) #folder10)
common_files = sorted(list(files1&files2&files3&files4&files5&files6&files7))

for filename in common_files:
    print(f"Comparing: {filename}")
    path1 = os.path.join(folder1,filename)
    path2 = os.path.join(folder2,filename)
    path3 = os.path.join(folder3,filename)
    path4 = os.path.join(folder4,filename)
    path5 = os.path.join(folder5,filename)
    path6 = os.path.join(folder11,filename)
    path7 = os.path.join(folder10,filename)
    
    try:
        X1 = process_data(path1)
        X2 = process_data(path2)
        X3 = process_data(path3)
        X4 = process_data(path4)
        X5 = process_data(path5)
        X6 = process_data(path6)
        X7 = process_data(path7)
        comparison_plot(X1,X2,X3,X4,X5,X6,X7, filename, output_dir1)
    except Exception as e:
        print(f"Error processing {filename}: {e}")
    
    
"""
    #----------------------------------------------------------------------
    #plot histograms plots
    #-----------------------------------------------------------------------
    
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
    #ax3.set_aspect('equal', adjustable='box')
    
    
    #----------------------------------------------------------------------
    #plot histograms plots
    #-----------------------------------------------------------------------
    
    
    fig2 = plt.figure(figsize=(12, 12))
    gs1 = gridspec.GridSpec(5, 5)
    
    # x 1D Histogram
    ax_x_hist = fig2.add_subplot(gs1[0, 0])
    ax_x_hist.hist(x1, bins=128, color='grey', histtype='step')
    ax_x_hist.set_xlim(-45, 45)
    ax_x_hist.set_xlabel("x [mm]")
    ax_x_hist.set_title("1D Hist x")
    ax_x_hist.spines["top"].set_visible(False)
    ax_x_hist.spines["right"].set_visible(False)
    
    # x-x'
    ax_x_xp = fig2.add_subplot(gs1[1, 0])
    ax_x_xp.hist2d(x1, xp1, bins=128, cmap='Greys')
    ax_x_xp.set_xlim(-45, 45)
    ax_x_xp.set_ylim(-4, 4)
    ax_x_xp.set_xlabel("x [mm]")
    ax_x_xp.set_ylabel("x' [mrad]")
    ax_x_xp.set_title("x' vs x")
    
    # x' 1D Histogram
    ax_xp_hist = fig2.add_subplot(gs1[1, 1])
    ax_xp_hist.hist(xp1, bins=128, color='grey', histtype='step')
    ax_xp_hist.set_xlim(-4, 4)
    ax_xp_hist.set_xlabel("x' [mrad]")
    ax_xp_hist.set_title("1D Hist x'")
    ax_xp_hist.spines["top"].set_visible(False)
    ax_xp_hist.spines["right"].set_visible(False)
    
    # x-y
    ax_x_y = fig2.add_subplot(gs1[2, 0])
    ax_x_y.hist2d(x1, y1, bins=128, cmap='Greys')
    ax_x_y.set_xlim(-45, 45)
    ax_x_y.set_ylim(-45, 45)
    ax_x_y.set_xlabel("x [mm]")
    ax_x_y.set_ylabel("y [mm]")
    ax_x_y.set_title("y vs x")
    
    # x'-y
    ax_xp_y = fig2.add_subplot(gs1[2, 1])
    ax_xp_y.hist2d(xp1, y1, bins=128, cmap='Greys')
    ax_xp_y.set_xlim(-4, 4)
    ax_xp_y.set_ylim(-45, 45)
    ax_xp_y.set_xlabel("x' [mrad]")
    ax_xp_y.set_ylabel("y [mm]")
    ax_xp_y.set_title("y vs x'")
    
    # y 1D Histogram
    ax_y_hist = fig2.add_subplot(gs1[2, 2])
    ax_y_hist.hist(y1, bins=128, color='grey', histtype='step')
    ax_y_hist.set_xlim(-45, 45)
    ax_y_hist.set_xlabel("y [mm]")
    ax_y_hist.set_title("1D Hist y")
    ax_y_hist.spines["top"].set_visible(False)
    ax_y_hist.spines["right"].set_visible(False)
    
    # y' vs x
    ax_yp_x = fig2.add_subplot(gs1[3, 0])
    ax_yp_x.hist2d(x1, yp1, bins=128, cmap='Greys')
    ax_yp_x.set_xlim(-45, 45)
    ax_yp_x.set_ylim(-4, 4)
    ax_yp_x.set_xlabel("x [mm]")
    ax_yp_x.set_ylabel("y' [mrad]")
    ax_yp_x.set_title("y' vs x")
    
    # y' vs x'
    ax_yp_xp = fig2.add_subplot(gs1[3, 1])
    ax_yp_xp.hist2d(xp1, yp1, bins=128, cmap='Greys')
    ax_yp_xp.set_xlim(-4, 4)
    ax_yp_xp.set_ylim(-4, 4)
    ax_yp_xp.set_xlabel("x' [mrad]")
    ax_yp_xp.set_ylabel("y' [mrad]")
    ax_yp_xp.set_title("y' vs x'")
    
    # y' vs y
    ax_yp_y = fig2.add_subplot(gs1[3, 2])
    ax_yp_y.hist2d(y1, yp1, bins=128, cmap='Greys')
    ax_yp_y.set_xlim(-45, 45)
    ax_yp_y.set_ylim(-4, 4)
    ax_yp_y.set_xlabel("y [mm]")
    ax_yp_y.set_ylabel("y' [mrad]")
    ax_yp_y.set_title("y' vs y")
    
    # y' 1D Histogram
    ax_yp_hist = fig2.add_subplot(gs1[3, 3])
    ax_yp_hist.hist(yp1, bins=128, color='grey', histtype='step')  # FIX: was y1 before
    ax_yp_hist.set_xlim(-4, 4)
    ax_yp_hist.set_xlabel("y' [mrad]")
    ax_yp_hist.set_title("1D Hist y'")
    ax_yp_hist.spines["top"].set_visible(False)
    ax_yp_hist.spines["right"].set_visible(False)
    
    # Save the figure
    plt.suptitle(f'Final Run Phase Space of {filename}', fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # reserve space for suptitle
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'Final_Run_1_Phase_Space_{filename}.png'))
    plt.close() 
"""

"""
    "ch10": cH10,
    "cV10": cV10,
    "cH12": cH12,
    "cV12": cV12,
    "cH13": cH13,
    "cV13": cV13,
    "cH14": cH14,
    "cV14": cV14,"""
        
        
    

