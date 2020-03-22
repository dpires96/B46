# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 03:28:45 2020

@author: Mohammad, Keval
"""

from sympy import *

def eigenvalues(rho, mass, V0):
    ###################################################
    ###################################################
    ########## LAMBDA SETUP ##########
    eig = symbols('eig')
    
    ##### Parameter Values #####
    S = 30.00 #[m^2]
    cbar = 2.0569 #[m]
    b = 15.911 #[m]
    
    MUc = mass/(rho*S*cbar)
    MUb = mass/(rho*S*b)
    
    Cx0, Cz0 = 1,1
    Cxu, Cxa, Cxq = -0.0279, -0.4797, -0.2817
    Czu, Cza, Czadot, Czq = -0.3762, -5.7434, -0.0035, -5.6629
    Cmu, Cma, Cmadot, Cmq, Ky = 0.0699, -0.5626, 0.178, -8.7941, sqrt(1.3925)
    
    CL = 2*mass/(rho*V0^2*S)
    CYB, CYBdot, CYp, CYr = -0.75, 0, -0.0304, 0.8495
    ClB, Clp, Kx, Clr, Kxz = -0.1026, -0.7108, sqrt(0.019), 0.2376, sqrt(0.002)
    CnB, CnBdot, Cnp, Cnr, Kz = 0.1348, 0, -0.0602, 0.2061, sqrt(0.042)
    ###################################################
    ###################################################
    ########## SYMMETRIC EOM :: EIGENVALUES ##########
    ##### Short Period Motion (2x2 Matrix) #####
    A_sSP = Matrix([
        [Cza + (Czadot-2*MUc)*eig, Czq + 2*MUc],
        [Cma + Cmadot*eig, Cmq - 2*MUc*Ky**2*eig]
        ])
    eig_sSP = solve(A_sSP.det(), eig)
    ##### Phugoid / Long Period Motion (3x3 Matrix) #####
    A_sLP = Matrix([
        [Cxu - 2*MUc*eig, Cz0, 0],
        [Czu, 0, MUc],
        [0, -eig, 1]
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
    
    ###################################################
    ###################################################
    ########## EXACT EIGENVALUES ##########
   
    
    print('eig_sSP', 'eig_sLP', 'eig_aAR', 'eig_aDR', 'eig_aAS')
    print(eig_sSP, eig_sLP, eig_aAR, eig_aDR, eig_aAS)
    return eig_sSP, eig_sLP, eig_aAR, eig_aDR, eig_aAS