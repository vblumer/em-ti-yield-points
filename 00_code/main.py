# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 10:12:13 2023

@author: BlumerVM
"""

## ----------------------------------------- DESCRIPTION -------------------------------------- ##
## Assemble and run a DAMASK simulation. Grain IDs and Euler angles are read from file created in
## Matlab MTEX. Material data is read from yaml file. Numerics settings and applied loading is 
## specified in the constants section.
## ----------------------------------------- END DESCRIPTION ---------------------------------- ##

import os
import numpy as np
import damask
import datetime
import pickle

## ----------------------------------------- UNITS     ---------------------------------------- ##
## ton, mm, s, N, MPa, Nmm 
## ----------------------------------------- END UNITS ---------------------------------------- ##

## ----------------------------------------- CONSTANTS ---------------------------------------- ##
FFT_iter            = 15        # Number of allowed FFT iterations. Standard is 250
num_cutbacks        = 0         # Number of allowed increment cutbacks. Standard is 3
num_of_yield_points = 9         # Number of stress paths / yield points
alpha_0             = 0         # First angle of stress decomposition in stress space (radians)
alpha_1             = np.pi     # Last angle of stress decomposition in stress space (radians)
applied_stress      = 2000      # Magnitude of the stress to be decomposed
N                   = 100       # Number of increments
t                   = 1         # length of time increment
## ----------------------------------------- END CONSTANTS ------------------------------------ ##

## ----------------------------------------- GRID --------------------------------------------- ##
## Read grain IDs from files created in MTEX, save grid as damask object
grid_pth                    = '../01_geom/grid_grain_IDs.txt'
dimensions_pth              = '../01_geom/dimensions.txt'

grid                        = np.loadtxt(grid_pth,delimiter=',').astype(int)
grid                        = grid-1                     # required because of mismatch in MTEX and DAMASK counting convention
grid                        = grid[..., np.newaxis]      #  M x N -> M x N x 1 dimension
cells                       = np.array(np.shape(grid))

[xmin,xmax,ymin,ymax,dx,dy] = np.loadtxt(dimensions_pth,delimiter=',')
physical_dimensions         = 1e-3*np.array([xmax-xmin,ymax-ymin,dx]) 
g                           = damask.Grid(grid, physical_dimensions)
## ----------------------------------------- END GRID ----------------------------------------- ##


## ----------------------------------------- MATERIAL------------------------------------------ ##
## Read Euler angles and material data from files, create damask material configuration object
euler_pth                                   = '../01_geom/grain_euler_angles.txt'
config_material                             = damask.ConfigMaterial()
config_material['homogenization']['dummy']  = {'N_constituents':1,'mechanical':{'type':'pass'}}
config_material['phase']['A']               = damask.ConfigMaterial.load('../02_mat/phenopower_somlo_frodal.yaml')

eulers                                      = np.loadtxt(euler_pth,delimiter=',')
O_A                                         = damask.Rotation.from_Euler_angles(eulers,degrees=True)

config_material                             = config_material.material_add(homogenization='dummy',phase='A',O=O_A)
## ----------------------------------------- END MATERIAL ------------------------------------- ##

## ----------------------------------------- WRITE SUBDIRECTORY ------------------------------- ##
## Create subdirectory where the analysis is run, change directory, save grid and material files
subdirectory_name   = 'results_{date:%Y-%m-%d_%H-%M-%S}'.format( date=datetime.datetime.now() )
directory_name      = '../99_results/'+subdirectory_name

os.makedirs(directory_name)
os.chdir(directory_name)

g.save('GRID.vti')
config_material.save('material.yaml')
## ----------------------------------------- END SUBDIRECTORY --------------------------------- ##

## ----------------------------------------- NUMERICS CONFIG ---------------------------------- ##
## Create file for damask numerics settings
with open('numerics.yaml', 'w') as f:
    f.write('grid:\n')
    f.write('  itmax: '     +str(FFT_iter)    +'\n')
    f.write('  maxCutBack: '+str(num_cutbacks)+'\n')
## ----------------------------------------- END NUMERICS ------------------------------------- ##

## ----------------------------------------- LOADCASE ----------------------------------------- ##
## Decompose stress, create damask loadcase,loadstep objects, save loadcase file
beta_space    = np.linspace(alpha_0,alpha_1,num_of_yield_points)

def decompose_xy(value, alpha): return [float(np.round(value*np.cos(alpha), 4)),
                                        float(np.round(value*np.sin(alpha),4))]

for i,beta in enumerate(beta_space):
    load_case         = damask.Config(solver={'mechanical':'spectral_basic'},loadstep=[])
    [sigma_1,sigma_2] = decompose_xy(applied_stress,beta)
    
    F           = [['x',   0,   0], 
                   ['x', 'x',   0], 
                   ['x', 'x', 'x']]
    P           = [[sigma_1,    'x', 'x'],
                   [     0,sigma_2, 'x'],
                   [     0,      0,   0]]
    
    loadstep    = {'boundary_conditions':{'mechanical':{'F':F,'P':P}},
                   'discretization':{'t':t,'N':N},'f_out':1}
    
    load_case['loadstep'].append(loadstep)
    
    load_case.save('LOADCASE_'+str(i)+'.yaml')
## ----------------------------------------- END LOADCASE ------------------------------------- ##

## ----------------------------------------- RUN JOB ------------------------------------------ ##
## Run jobs for all stress paths, record runtimes
for idx,beta in enumerate(beta_space):    
    start   = datetime.datetime.now()
    
    os.system('DAMASK_grid --geom GRID.vti --load LOADCASE_'+str(idx)+'.yaml')
    
    runtime = datetime.datetime.now()-start
    with open('runtime_'+str(idx)+'.txt', 'w') as text_file:
        text_file.write(str(runtime.seconds)+'s')
                                            

    ## -- POSTPROCESS -- ##
    ## Read hdf5, add Green Lagrange strain E, 2nd Piola Kirchhoff stress S
    r       = damask.Result('GRID_LOADCASE_'+str(idx)+'.hdf5')
    
    E       = r.get('epsilon_V^1(F)')
    if E is None:
        r.add_strain(m=1)
        E = r.get('epsilon_V^1(F)')
    
    S = r.get('S')
    if S is None:
        r.add_stress_second_Piola_Kirchhoff()
        S = r.get('S')
    
    ## Obtain increment identifiers, initialize arrays
    keys       = list(E.keys())
    
    E_macro    = np.zeros([3,3,len(keys)])
    S_macro    = np.zeros([3,3,len(keys)])
    
    E_eq_macro = np.zeros(len(keys))
    S_eq_macro = np.zeros(len(keys))
    
    ## Compute volume averaged E,S and equivalent sqrt(E:E), sqrt(S:S),
    for i,key in enumerate(keys):
        E_macro[:,:,i] = np.mean(E[key],axis=0)
        S_macro[:,:,i] = np.mean(S[key],axis=0)
        
        E_eq_macro[i]  = np.sqrt(np.tensordot(E_macro[:,:,i],E_macro[:,:,i]))
        S_eq_macro[i]  = np.sqrt(np.tensordot(S_macro[:,:,i],S_macro[:,:,i]))
    
    ## Save relevant data to pickle, delete hdf5 file 
    with open('results_'+str(idx)+'.pickle', 'wb') as handle:
        pickle.dump([E_macro,S_macro,E_eq_macro,S_eq_macro], 
                    handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    os.remove('GRID_LOADCASE_'+str(idx)+'.hdf5')
## ----------------------------------------- END RUN JOB -------------------------------------- ##
