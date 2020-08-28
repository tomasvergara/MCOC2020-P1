#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 23:09:31 2020

@author: Tomas
"""

import scipy as sp
from scipy.integrate import odeint
from numpy import matrix , cos , sin , transpose , array , zeros 
import matplotlib.pylab as plt


def satelite (z,t):
    R = array([[ cos(om*t) , -sin(om*t) , 0 ] ,
               [ sin(om*t) , cos(om*t ) , 0 ] ,
               [ 0 , 0 , 1]])

    Rt = transpose(R)

    R_1 = array([[ -sin(om*t)*om , -cos(om*t)*om , 0 ] ,
                 [ cos(om*t)*om , -sin(om*t )*om , 0 ] ,
                 [ 0 , 0 , 0]])

    R_2 = array([[ -cos(om*t)*om**2 , sin(om*t)*om**2 , 0 ] ,
                 [ -sin(om*t)*om**2 , -cos(om*t )*om**2 , 0 ] ,
                 [ 0 , 0 , 0]])
    
    zp = zeros(6)
    zp[0:3] = z[3:6]
    zp[3:6] = (-G * m_t * z[0:3])/(r**3) - Rt @  (R_2@z[0:3] + 2*R_1@z[3:6] )
    return zp

G = 6.67 * (10**-11) # m^3/(kg*s^2)
m_t = 5.972 * (10**24) # kg
d = 700000 # m
r = 7071000
atm = 80000

x = 6371000 + d # m
y = 0
z = 0
vx = 0
vy = 7000
vz = 0
z0 = array([ x , y , z , vx , vy , vz ])


pi = 3.1415

# Vector de tiempo
t = sp.linspace( 0 , 86400 , 1001)

om = 2 * pi / (24*3600)

sol = odeint( satelite , z0 , t)

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

plt.figure()
plt.grid()
plt.plot(t,x)
plt.plot(t,y)
plt.plot(t,z)

plt.title("Historias de tiempo")
plt.ylabel("Distancia en miles de km")
plt.xlabel("Tiempo [s]")
plt.legend(["X (t)", "Y (t)" , "Z (t)"])
plt.xlim(0,12700)
plt.savefig("grafico_vstiempo", bbox_inches = "tight")