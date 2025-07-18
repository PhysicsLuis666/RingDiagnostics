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

def density(file_path2):
    
    df = pd.read_csv(file_path2, sep="\s+",comment="%",header=None, skip_blank_lines=True)
    
    rho = df.to_numpy()
    
    return rho

def extract_number(filename):
    """Extract 4-digit number from filename, e.g., 0050 from 'bunch_0050.txt'."""
    match = re.search(r'(\d{4})', filename)
    return match.group(1) if match else None

def emittance(X1):
    x1, xp1, y1, yp1 = X1[:, 0], X1[:, 1], X1[:, 2], X1[:, 3] #FR1
    covx = np.cov(x1,xp1)
    covy = np.cov(y1,yp1)
    ex = np.sqrt(np.var(x1)*np.var(xp1)-covx[0,1]**2)
    ey = np.sqrt(np.var(y1)*np.var(yp1)-covy[0,1]**2)
    
    return ex,ey
    

def comparison_plot(X1,rho, filename, output_dir): 
    x1, xp1, y1, yp1 = X1[:, 0], X1[:, 1], X1[:, 2], X1[:, 3] #FR1
    nx = np.sum(np.array(rho), axis=1)
    ny = np.sum(np.array(rho), axis=0)
    
    nx_norm = nx/np.sum(nx)
  
    
    beamx = np.std(x1)
    beamy = np.std(y1)
    
    rx = 2*beamx
    ry = 2*beamy
    ns = len(x1)/(np.pi*rx*ry)
    print(f"ns = {ns}")
    print(f"rx = {rx} , ry = {ry}")

    Nx = 2*ns*np.sqrt(1-(x1/rx)**2)
    Ny = 2*ns*np.sqrt(1-(y1/ry)**2)
    
    xmax = np.max(x1)
    xmin = np.min(x1)
    ymax = np.max(y1)
    ymin = np.min(y1)
    print(f"max x: {xmax}")
    print(f"max y: {ymax}")
    
    
    dx = (xmax-xmin)/(128-1)
    xcoord = xmin + dx * np.arange(128)
    xcoord -= np.mean(xcoord)
    
    histogram_x1, edges = np.histogram(x1,  bins=128, range=(xcoord[0], xcoord[-1]))
    x_bins = 0.5 * (edges[:-1] + edges[1:])

    
    dy = (ymax-ymin)/(128-1)
    ycoord = ymin + dy*np.arange(128)
    ycoord -= np.mean(ycoord)
    
  
    fig1, (ax1,ax2, ax3) = plt.subplots(3,1, figsize = (6,12))
    ax1.plot(xcoord, nx, label='From 2D rho grid',  linestyle="--")
   #ax1.set_xlim(-xmax, xmax)
    #ax1.set_ylim(0, 3)
    ax1.set_xlabel("x [mm]")
    ax1.set_ylabel("n(x)")
    ax1.grid(True)
    ax1.legend()
    #ax1.set_ylim(ymax,ymin)
    ax2.plot(ycoord, ny)
   # ax2.set_xlim(-ymax,ymax)
    #ax2.set_ylim(0, 3)
    ax2.set_xlabel("y [mm]")
    ax2.set_ylabel("n(y)")
    ax2.grid(True)
    
    print("area:", np.sum(nx)*dx)
    #ax3.hist(x1,bins=128, label="From Raw Data", histtype="step", color="green")
    
    #ax3.plot(xcoord, nx_norm, label = "normalized nx")
    #ax3.grid(True)
    #ax1.set_ylim(0, 3)
    #ax3.legend()
    
    ax3.hist(x1, bins=128, label="raw data histogram", histtype="step", color="green")
    ax3.set_ylabel("counts")
    ax3.set_xlabel("x [mm]")
    ax3.grid(True)
    ax3.legend()
    
    save_path = os.path.join(output_dir1, f"histogram_{filename}.png")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
    
    centroid_x = np.sum(xcoord * nx) / np.sum(nx)
    print(f"Centroid of x density: {centroid_x:.6e} m")
    

    
  
        
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
    """
    plt.suptitle(f"1D Histogram Comparisons at Different Kicker Settings")
    plt.tight_layout()
    os.makedirs(output_dir3, exist_ok=True)
    plt.savefig(os.path.join(output_dir3, f'ComparingFilesLargeSIVpoint.99toLowerSIVpoints_{filename}.png'))
    plt.close(fig4)
    """
    
    
        
        
    
      

folder1 = "outputs3-1/TrialRun2"  #folder run .418 sfv and .77 siv
folder2 = "outputs3-2/density2"
output_dir1 = "/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs3-1/Analysis1-3"
#output_dir2 = "/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs2/Analysis1-4"
#output_dir3 = "/Users/luisruiz/sns-ring-model/scripts/full_injection_benchmark/outputs2/Analysis1-5"


# Get files and map them by the extracted number
files1 = os.listdir(folder1)
files2 = os.listdir(folder2)

filemap1 = {extract_number(f): f for f in files1 if extract_number(f)}
filemap2 = {extract_number(f): f for f in files2 if extract_number(f)}

# Find common numbers
common_numbers = sorted(set(filemap1) & set(filemap2))

for number in common_numbers:
    f1 = filemap1[number]
    f2 = filemap2[number]

    path1 = os.path.join(folder1, f1)
    path2 = os.path.join(folder2, f2)

    print(f"Comparing files for turn {number}: {f1} & {f2}")

    try:
        X1 = process_data(path1)
        Rho = density(path2)
        comparison_plot(X1, Rho, f"turn_{number}", output_dir1)
    except Exception as e:
        print(f"Error processing turn {number}: {e}")
    
    
