# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 03:28:45 2020

@author: Mohammad, Keval
"""

from sympy import *

def eigenvalues(mass, V0):
    ###################################################
    ###################################################
    ########## LAMBDA SETUP ##########
    eig = symbols('eig')
    
    ##### Parameter Values #####
    S = 24.2 #[m^2]
    cbar = 2.022 #[m]
    b = 13.36 #[m]
    
    MUc = 102.7
    MUb = 15.5
    
    Cx0, Cz0 = 0, -1.136
    Cxu, Cxa, Cxq = -0.2199, 0.4653, 0
    Czu, Cza, Czadot, Czq = -2.272, -5.16, -1.43, -3.86
    Cmu, Cma, Cmadot, Cmq, Ky = 0, -0.43, -3.7, -7.04, sqrt(0.98)
    
    CL = 1.136
    CYB, CYBdot, CYp, CYr = -0.9896, 0, -0.087, 0.43
    ClB, Clp, Kx, Clr, Kxz = -0.0772, -0.3444, sqrt(0.012), 0.28, 0.002
    CnB, CnBdot, Cnp, Cnr, Kz = 0.1638, 0, -0.0108, -0.193, sqrt(0.037)
    ###################################################
    ###################################################
    ########## SYMMETRIC EOM :: EIGENVALUES ##########
    ##### Short Period Motion (2x2 Matrix) #####
    A_sSP = Matrix([
        [Cza + (Czadot-2*MUc)*eig, Czq + 2*MUc],
        [Cma + Cmadot*eig, Cmq - 2*MUc*Ky**2*eig]
        ])
    eig_sSP = solve(A_sSP.det(), eig)
    ##### Phugoid / Long Period Motion (4x4 Matrix) #####
    A_sLP = Matrix([
        [Cxu - 2*MUc*eig, Cxa, Cz0, 0],
        [Czu, Cza, 0, Czq + 2*MUc],
        [0, 0, -eig, 1],
        [Cmu, Cma, 0, Cmq]
        ])
    eig_sLP = solve(A_sLP.det(), eig)
    
    ###################################################
    ###################################################
    ########## ASYMMETRIC EOM :: EIGENVALUES ##########
    ##### Aperiodic Roll Motion (1x1 Matrix) #####
    A_aAR = Matrix([
        [Clp - 4*MUb*Kx**2*eig]
        ])
    eig_aAR = solve(A_aAR.det(), eig)
    ##### Dutch Roll Motion (2x2 Matrix) #####
    A_aDR = Matrix([
        [CYB - 2*MUb*eig, -4*MUb],
        [CnB, Cnr - 4*MUb*Kz**2*eig]
        ])
    eig_aDR = solve(A_aDR.det(), eig)
    ##### Aperiodic Spiral Motion (1x1 Matrix) #####
    A_aAS = Matrix([
        [CYB, CL, 0, -4*MUb],
        [0, -0.5*eig, 1, 0],
        [ClB, 0, Clp, Clr],
        [CnB, 0, Cnp, Cnr]
        ])
    eig_aAS = solve(A_aAS.det(), eig)
    ##### Dutch Roll & Aperiodic Roll Motion (1x1 Matrix) #####
    A_aDRAR = Matrix([
        [-ClB + 0.5*Clr*eig + 2*MUb*Kxz*eig**2, Clp - 4*MUb*Kx**2*eig],
        [-CnB + 0.5*Cnr*eig - 2*MUb*Kz**2*eig**2, Cnp + 4*MUb*Kxz*eig]
        ])
    eig_aDRAR = solve(A_aDRAR.det(), eig)
    
    ###################################################
    ###################################################
    ########## EXACT EIGENVALUES ##########
   
    
    print('eig_sSP', 'eig_sLP', 'eig_aAR', 'eig_aAS', 'eig_aDR')
    print(eig_sSP, eig_sLP, eig_aAR, eig_aAS, eig_aDR)
    return eig_sSP, eig_sLP, eig_aAR, eig_aAS, eig_aDR



mass = 4547.8
V0 = 59.9

eigenvalues(mass, V0)
