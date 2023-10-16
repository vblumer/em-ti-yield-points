#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 11:29:41 2023

@author: blumervm
"""

## ----------------------------------------- DESCRIPTION -------------------------------------- ##
## Postprocess DAMASK analyses. From the specified folder, all recorded pickle files are read.
## For each stress path, the yield point is identified using the 0.2% offset method. 
## ----------------------------------------- END DESCRIPTION ---------------------------------- ##

import numpy as np 
import pickle as pkl
import glob
import matplotlib.pyplot as plt

## ----------------------------------------- READ --------------------------------------------- ##
pth = '../99_results/results_2023-09-22_16-23-27/'

loadcase_files = glob.glob(pth+'*.pickle')
loadcase_files = [w.split('.pickle')[0].split('_')[-1] for w in loadcase_files]
## ----------------------------------------- END READ ----------------------------------------- ##

## ----------------------------------------- COMPUTE ------------------------------------------ ##
index          = np.sort(np.array([loadcase_files][0]).astype(int))
yield_points   = np.zeros([len(index),2])
for j,lc in enumerate(index):
    # Read results of one stress path
    path = pth+'results_'+str(j)+'.pickle'
    with open(path, 'rb') as handle:
        [E_macro,S_macro,E_eq_macro,S_eq_macro] = pkl.load(handle)
    
    # Compute index of intersection of stress/strain curve and linear offset
    lin_stiffness     = (S_eq_macro[1]-S_eq_macro[0])/(E_eq_macro[1]-E_eq_macro[0])
    E_eq_macro_offset = E_eq_macro-0.002
    y                 = lin_stiffness*E_eq_macro_offset
    idx               = np.where(S_eq_macro<y)[0]
    
    # If length of idx is zero, than the criterion is not fulfilled at any stage
    if len(idx)!=0:
        idx       = idx[0]
        # Take average of values before and after intersection
        yld_pt_eq = [np.mean([E_eq_macro[idx],E_eq_macro[idx-1]]),
                     np.mean([S_eq_macro[idx],S_eq_macro[idx-1]])]
        
        # Plot stress/strain curve
        plt.figure(dpi=1200)
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        
        ax = plt.gca()
        ax.grid()
        ax.set_axisbelow(True)
        plt.plot(E_eq_macro, S_eq_macro)
        plt.plot(E_eq_macro, y)
        plt.scatter(yld_pt_eq[0],yld_pt_eq[1])
        plt.show()
        
        # Record yield point
        yield_points[j,:] = np.array([np.mean([S_macro[0,0,idx],S_macro[0,0,idx-1]]),
                                      np.mean([S_macro[1,1,idx],S_macro[1,1,idx-1]])])
    else:
        print('No yielding detected')
## ----------------------------------------- END COMPUTE -------------------------------------- ##

## ----------------------------------------- PLOT --------------------------------------------- ##
plt.figure(dpi=1200)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
xmax =  1400
xmin = -1400
plt.xticks(np.linspace(xmin,xmax,15))

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.grid()
ax.set_axisbelow(True)

for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)

plt.scatter(1e-6*yield_points[:,0],1e-6*yield_points[:,1],c='#ec7900')

plt.savefig('../98_pictures/plot_yld.png', bbox_inches='tight')
plt.show()
## ----------------------------------------- END PLOT ----------------------------------------- ##
