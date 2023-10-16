#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 16:46:53 2023

@author: blumervm
"""

import numpy as np
import damask
import matplotlib.pyplot as plt
import glob
pth = '../99_results/results_2023-09-15_11-41-33/'

loadcase_files = glob.glob(pth+'LOADCASE_*')
loadcase_files = [w.split('LOADCASE_')[-1].split('.')[0] for w in loadcase_files]

index          = np.sort(np.array([loadcase_files][0]).astype(int))

yield_points = np.zeros([len(index),2])
for j,lc in enumerate(index):
    path = pth+'GRID_LOADCASE_'+str(lc)+'.hdf5'
    
    r = damask.Result(path)
    
    E = r.get('epsilon_V^1(F)')
    if E is None:
        r.add_strain(m=1)
        E = r.get('epsilon_V^1(F)')
    
    S = r.get('S')
    if S is None:
        r.add_stress_second_Piola_Kirchhoff()
        S = r.get('S')
    
    keys = list(E.keys())
    
    E_macro = np.zeros([3,3,len(keys)])
    S_macro = np.zeros([3,3,len(keys)])
    
    E_eq_macro = np.zeros(len(keys))
    S_eq_macro = np.zeros(len(keys))
    for i,key in enumerate(keys):
        E_macro[:,:,i] = np.mean(E[key],axis=0)
        S_macro[:,:,i] = np.mean(S[key],axis=0)
        
        E_eq_macro[i] = np.sqrt(np.tensordot(E_macro[:,:,i],E_macro[:,:,i]))
        S_eq_macro[i] = np.sqrt(np.tensordot(S_macro[:,:,i],S_macro[:,:,i]))
    
    
    lin_stiffness = (S_eq_macro[1]-S_eq_macro[0])/(E_eq_macro[1]-E_eq_macro[0])
    E_eq_macro_offset = E_eq_macro-0.002
    y = lin_stiffness*E_eq_macro_offset
    
    idx    = np.where(S_eq_macro<y)[0]
    if len(idx)!=0:
        idx = idx[0]
        yld_pt_eq =  [np.mean([E_eq_macro[idx],E_eq_macro[idx-1]]),
                      np.mean([S_eq_macro[idx],S_eq_macro[idx-1]])]
        
        plt.plot(E_eq_macro, S_eq_macro)
        plt.plot(E_eq_macro, y)
        plt.scatter(yld_pt_eq[0],yld_pt_eq[1])
        plt.show()
        
        yield_points[j,:] = np.array([np.mean([S_macro[0,0,idx],S_macro[0,0,idx-1]]),
                                      np.mean([S_macro[1,1,idx],S_macro[1,1,idx-1]])])
    else:
        print('No yielding detected')

plt.scatter(yield_points[:,0],yield_points[:,1])
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.show()
