# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 03:28:45 2020

@author: Mohammad, Keval
"""

from sympy import *
from math import sin, cos

###################################################
###################################################
########## LAMBDA SETUP ##########
eig = symbols('eig')

##### Parameter Values #####
# S = 24.2 #[m^2]
# cbar = 2.022 #[m]
# b = 13.36 #[m]
# V = 59.9 #[m/s]

# MUc = 102.7
# MUb = 15.5

# Cx0, Cz0 = 0, -1.136
# Cxu, Cxa, Cxq = -0.2199, 0.4653, 0
# Czu, Cza, Czadot, Czq = -2.272, -5.16, -1.43, -3.86
# Cmu, Cma, Cmadot, Cmq, Ky = 0, -0.43, -3.7, -7.04, sqrt(0.98)

# CL = 1.136
# CYB, CYBdot, CYp, CYr = -0.9896, 0, -0.087, 0.43
# ClB, Clp, Kx, Clr, Kxz = -0.0772, -0.3444, sqrt(0.012), 0.28, 0.002
# CnB, CnBdot, Cnp, Cnr, Kz = 0.1638, 0, -0.0108, -0.193, sqrt(0.037)
# ##### Parameter Values #####
S = 30.00 #[m^2]
cbar = 2.0569 #[m]
b = 15.911 #[m]

V = [86.90426712,
90.27120978,
85.80364147,
93.55135357,
88.83844801,
90.17191004]
rho = [1.033917722,
1.035874106,
1.034387703,
1.036636471,
1.035305677,
0.972332419]
mass = [6266.730308,
6229.230633,
6218.934094,
6209.091148,
6203.058374,
6177.067553]
th0 = [0.147148709,
0.054197709,
0.090490086,
0.060699061,
0.05714208,
0.085618872]

i = 2

MUc = mass[i]/(rho[i]*S*cbar)
MUb = mass[i]/(rho[i]*S*b)
CL = 2*mass[i]*9.81/(rho[i]*V[i]**2*S)
Cx0 = 2*mass[i]*9.81*sin(th0[i])/(rho[i]*V[i]**2*S)
Cz0 = -2*mass[i]*9.81*cos(th0[i])/(rho[i]*V[i]**2*S)

Cxu, Cxa, Cxq = -0.0279, -0.4797, -0.2817
Czu, Cza, Czadot, Czq = -0.3762, -5.7434, -0.0035, -5.6629
Cmu, Cma, Cmadot, Cmq, Ky = 0.0699, -0.5626, 0.178, -8.7941, sqrt(1.3925)

CYB, CYBdot, CYp, CYr = -0.75, 0, -0.0304, 0.8495
ClB, Clp, Kx, Clr, Kxz = -0.1026, -0.7108, sqrt(0.019), 0.2376, 0.002
CnB, CnBdot, Cnp, Cnr, Kz = 0.1348, 0, -0.0602, -0.2061, sqrt(0.042)
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
i = 0

MUc = mass[i]/(rho[i]*S*cbar)
MUb = mass[i]/(rho[i]*S*b)
CL = 2*mass[i]*9.81/(rho[i]*V[i]**2*S)
Cx0 = 2*mass[i]*9.81*sin(th0[i])/(rho[i]*V[i]**2*S)
Cz0 = -2*mass[i]*9.81*cos(th0[i])/(rho[i]*V[i]**2*S)
A_sLP = Matrix([
    [Cxu - 2*MUc*eig, Cxa, Cz0, 0],
    [Czu, Cza, 0, Czq + 2*MUc],
    [0, 0, -eig, 1],
    [Cmu, Cma, 0, Cmq]
    ])
# A_sLP = Matrix([
#     [Cxu - 2*MUc*eig, Cz0, 0],
#     [Czu, 0, 2*MUc],
#     [0, -eig, 1]
#     ])
eig_sLP = solve(A_sLP.det(), eig)

###################################################
###################################################
########## ASYMMETRIC EOM :: EIGENVALUES ##########
##### Aperiodic Roll Motion (1x1 Matrix) #####
i = 1

MUc = mass[i]/(rho[i]*S*cbar)
MUb = mass[i]/(rho[i]*S*b)
CL = 2*mass[i]*9.81/(rho[i]*V[i]**2*S)
Cx0 = 2*mass[i]*9.81*sin(th0[i])/(rho[i]*V[i]**2*S)
Cz0 = -2*mass[i]*9.81*cos(th0[i])/(rho[i]*V[i]**2*S)
A_aAR = Matrix([
    [Clp - 4*MUb*Kx**2*eig]
    ])
eig_aAR = solve(A_aAR.det(), eig)
##### Dutch Roll Motion (2x2 Matrix) #####
i = 3

MUc = mass[i]/(rho[i]*S*cbar)
MUb = mass[i]/(rho[i]*S*b)
CL = 2*mass[i]*9.81/(rho[i]*V[i]**2*S)
Cx0 = 2*mass[i]*9.81*sin(th0[i])/(rho[i]*V[i]**2*S)
Cz0 = -2*mass[i]*9.81*cos(th0[i])/(rho[i]*V[i]**2*S)
A_aDR = Matrix([
    [CYB - 2*MUb*eig, -4*MUb],
    [CnB, Cnr - 4*MUb*Kz**2*eig]
    ])
eig_aDR = solve(A_aDR.det(), eig)
##### Aperiodic Spiral Motion (4x4 Matrix) #####
i = 5

MUc = mass[i]/(rho[i]*S*cbar)
MUb = mass[i]/(rho[i]*S*b)
CL = 2*mass[i]*9.81/(rho[i]*V[i]**2*S)
Cx0 = 2*mass[i]*9.81*sin(th0[i])/(rho[i]*V[i]**2*S)
Cz0 = -2*mass[i]*9.81*cos(th0[i])/(rho[i]*V[i]**2*S)
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

#Validating analytical model
# E1 = -0.039161 + 0.037971*I
# E2 = -0.00029107 + 0.0066006*I
# E3 = -0.32910
# E4 = -0.01080
# E5 = -0.058028 + 0.37700*I

#Verifying numerical model
E1 = -0.033292 + 0.063521*I
E2 = -0.000066201 + 0.0034606*I
E3 = -0.81614
E4 = 0.0019494
E5 = -0.059847 + 0.38524*I

A = abs(E1)
B = abs(eig_sSP[0])
print((A-B)/max(A,B)*100)

A = abs(E2)
B = abs(eig_sLP[0])
print((A-B)/max(A,B)*100)

A = abs(E3)
B = abs(eig_aAR[0])
print((A-B)/max(A,B)*100)

A = abs(E4)
B = abs(eig_aAS[0])
print((A-B)/max(A,B)*100)

A = abs(E5)
B = abs(eig_aDR[0])
print((A-B)/max(A,B)*100)
