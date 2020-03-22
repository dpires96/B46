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


##### Short Period Motion #####
C1symsp = Matrix([
    [(Czadot-2*MUc)*cbar/V, 0],
    [Cmadot*cbar/V, -2*MUc*Ky**2*cbar/V]
    ])

C2symsp = Matrix([
    [Cza, Czq + 2*MUc],
    [Cma, Cmq]
    ])

Asymsp = -C1symsp**-1*C2symsp
Asymsp = Asymsp - eig*eye(max(Asymsp.shape))

solve(Asymsp.det(), eig)