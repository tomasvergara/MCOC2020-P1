#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:48:11 2020

@author: Tomas
"""

import numpy as np
from numpy import array , exp , zeros , cos , pi
from matplotlib.pylab import linspace
from scipy.integrate import odeint
import matplotlib.pylab as plt

# m * x'' + c * x' + k * x = 0  

# CONDICIONES INICIALES
x0 = 1
v0 = 1
z0 = [x0,v0]

# PARAMETROS
m = 1.0 # kg
f = 1.0 # Hz
e = 0.2
w = 2 * pi * f
k = m * (w**2)
c = 2 * e * w * m
t = linspace(0 , 4. , 100)

def z_punto(z , t):
    zp = np.zeros(2)
    zp[0] = z[1]
    z1 = z[0]
    z2 = z[1]
    zp[1] = -(c * z2 + k * z1)/m
    return zp

def zp_euler(zp , z0 , t , Nsubdivisiones):
    Nt = len(t)
    Ndim = len(array(z0))
    z = zeros((Nt,Ndim))
    z[0,:] = z0[0]
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsubdivisiones
        z_temp = z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp += dt * z_punto(z_temp,t_anterior + k*dt)
        z[i,:] = z_temp
    return z

# GRAFICO
analitica = exp(-c*t/2) * cos(w*t)
plt.plot(t , analitica , "k" , linewidth=2 ,label = "sol analitica")

sol = odeint(z_punto,z0,t)
             
x_2punto = sol[:,0]
plt.plot(t , x_2punto , "-b" , label = "odeint")

## Grafico con distintas subdivisiones
zp = sol[:,:]
sol1 = zp_euler(z_punto , z0 , t , 1)
euler1 = sol1[:,0]
plt.plot(t , euler1 , "--g" , label = "Nsub = 1")

sol2 = zp_euler(z_punto , z0 , t , 10)
euler2 = sol2[:,0]
plt.plot(t , euler2 , "--r" , label = "Nsub = 10")

sol3 = zp_euler(z_punto , z0 , t , 100)
euler3 = sol3[:,0]
plt.plot(t , euler3 , "--" , color = 'orange' , label = "Nsub = 100")

plt.title("Oscilador armónico")
plt.xlabel("Tiempo [s]")
plt.ylabel("X (t)")
plt.legend()
plt.axis([0,4,-1.2,1.2])
plt.savefig("oscilador_armonico.png" , bbox_inches = "tight")