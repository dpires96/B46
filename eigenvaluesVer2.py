# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 03:28:45 2020

@author: Mohammad, Keval
"""

from sympy import *


########## LAMBDA SETUP ##########
eig = symbols('eig')

########## SYMMETRIC EOM TO STATE SPACE ##########
##### Parameters #####
Cxu, MUc, Cxa, Cz0, Cxq = symbols('Cxu, MUc, Cxa, Cz0, Cxq')

Czu, Cza, Czadot, Cx0, Czq = symbols('Czu, Cza, Czadot, Cx0, Czq')

Cmu, Cma, Cmadot, Cmq, Ky = symbols('Cmu, Cma, Cmadot, Cmq, Ky')

cbar, V = symbols('cbar, V')


##### Short Period Motion (2x2 Matrix) #####
Asymsp = Matrix([
    [Cza + (Czadot-2*MUc)*eig, Czq + 2*MUc],
    [Cma + Cmadot*eig, Cmq - 2*MUc*Ky**2*eig]
    ])

eigsymsp = solve(Asymsp.det(), eig)