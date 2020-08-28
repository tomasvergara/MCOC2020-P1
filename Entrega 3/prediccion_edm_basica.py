#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:56:46 2020

@author: Tomas
"""

import scipy as sp
from scipy.integrate import odeint
from numpy import cos , sin , transpose , array , zeros 
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
r = 7071000


from datetime import datetime 

ti = "2020-07-22T22:59:42.000000"
ti = ti.split("T")
ti = "{} {}".format( ti[0] , ti[1] )
ti = datetime.strptime( ti, "%Y-%m-%d %H:%M:%S.%f" )

tf = "2020-07-24T00:59:42.000000"
tf = tf.split("T")
tf = "{} {}".format( tf[0] , tf[1] )
tf = datetime.strptime( tf , "%Y-%m-%d %H:%M:%S.%f" )

deltaT = ( tf - ti ).total_seconds()

xi = -2054124.071387
yi = 6466609.881937
zi = 2001866.224693
vxi = 893.626490
vyi = 2499.053195
vzi = -7119.273313

xf = -948622.852603
yf = 346929.024224
zf = 6992854.262127
vxf = 1957.902881
vyf = 7324.709077
vzf = -97.529820

z0 = array([ xi , yi , zi , vxi , vyi , vzi ])


pi = 3.1415

# Vector de tiempo
t = sp.linspace( 0 , deltaT , 9361)

om = 2 * pi / ((24)*3600)

sol = odeint( satelite , z0 , t)

x = sol[:,:]

pos_final = array([ xf , yf , zf , vxf , vyf , vzf ]) - sol[-1]

norma_distancia_error = (pos_final[0]**2 + pos_final[1]**2 + pos_final[2]**2)**(1/2.0)

print (norma_distancia_error)
    
    







