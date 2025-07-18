#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 11:03:44 2025

@author: luisruiz
"""
"""
test script for waveform

"""
import numpy as np
import matplotlib.pyplot as plt

def strength(ti, time, tf, si, sf):
    strength_array = np.empty_like(time)

    mask_before = time < ti
    strength_array[mask_before] = si

    mask_after = time > tf
    strength_array[mask_after] = sf

    mask_between = ~(mask_before | mask_after)
    dt = np.sqrt((time[mask_between] - ti) / (tf - ti))
    strength_array[mask_between] = (si - sf) * (1 - dt) + sf

    return strength_array

# Time array: from -1 ms to 1 ms

# Define waveform parameters
tih = -0.001                  # unshifted waveform
tiv = -.0004 #-0.000484             # shifted waveform starts earlier
tfinal = 0.001
sih = 0.99
sfh = 0.43
siv = 0.56
sfv = 0.12

if tih < 0.0 or tiv < 0.0:
    time = np.linspace(-0.001, 0.001, 1000)
elif tih >0.0 or tiv >0.0:
    time = np.linspace(0.0, 0.001, 1000)
else:
    time = np.linspace(0.0,.001,1000)
    

# Calculate waveforms
horizontalS = strength(tih, time, tfinal, sih, sfh)
verticalS   = strength(tiv, time, tfinal, siv, sfv)
verticalWaveform2 =  strength(0.000, time, tfinal, .99, sfv)
# Plot
plt.plot(time * 1000, verticalWaveform2, label=f"vertical Waveform 2 (tiv = .99)")
plt.plot(time * 1000, verticalS, label=f"Vertical Waveform 1(tiv ={tiv})")
plt.axvline(0.0, color="gray", linestyle="--")
plt.xlabel("Time [ms]")
plt.ylabel("Strength")
plt.title("Time-Shifted Square Root Waveforms")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()