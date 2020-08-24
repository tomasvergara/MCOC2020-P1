#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 09:44:49 2020

@author: Tomas
"""

import scipy as sp
from scipy.integrate import odeint

# Unidades base

cm = 0.01
inch = 2.54 * cm
g = 9.81
m = 15. #kg

# Viento
V = 20. #m/s


# Coeficiemte de arrastre

rho = 1.225 #kg / m**3
cd = 0.47
D = 8.5 * inch
r = D/2
A = sp.pi * r**2
CD = 0.5 * rho * cd * A 

# Funcion a intregrar
# z vector de estado = [x , y, vx, vy]
# dz/dt = bala(z,t)
#       [    z2     ]
# dz/dt=[           ]
#       [  FD/m  - g  ]

# Vector de estado
# z[0]   -> x
# z[1]   -> y
# z[2]   -> vx
# z[3]   -> vy

def bala(z,t):
    zp = sp.zeros(4)
    zp[0] = z[2]
    zp[1] = z[3]
    v = z[2:4]
    v[0] = v[0] - V
    vnorm = sp.sqrt(sp.dot(v,v))
    FD = ( -CD * sp.dot(v,v) * v ) / (vnorm)
    zp[2] = FD[0]/ m
    zp[3] = ( FD[1]/ m ) - g
    
    return zp


# Vector de tiempo
t = sp.linspace(0 , 5.6 , 1001)

# Parte en el origen 
vi = 100*1000/3600.
z0 = sp.array([0,0,vi,vi])

sol = odeint( bala , z0 , t)

import matplotlib.pylab as plt

x = sol[:,0]
y = sol[:,1]

plt.figure(1)
plt.plot(x,y)
plt.axis([0,50,0,50])
plt.show()
