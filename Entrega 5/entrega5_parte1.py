#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:54:22 2020

@author: Tomas
"""

import  matplotlib.pylab as plt
import numpy as np
from scipy import pi , sin , cos , dot , sqrt , zeros
from scipy.integrate import odeint
from leereof import leer_eof

# PARAMETROS
G = 6.67 * (10**-11) 
m_t = 5.98 * (10**24) 
omega = (2*pi) / 86400 

# FUNCION
def zpunto(z,t):
    zp = zeros(6)
    sin_o = sin(omega*t)
    cos_o = cos(omega*t)
    z1 = z[0:3]
    r = sqrt(dot(z1,z1))
    
    RT = np.array([[cos_o,sin_o,0],
                  [-sin_o,cos_o,0],
                  [0,0,1]])
    
    R1 = np.array([[-sin_o*omega,-cos_o*omega,0],
                   [cos_o*omega, -sin_o*omega,0],
                   [0,0,0]])  
    
    R2 = np.array([[-cos_o*(omega**2),sin_o*(omega**2),0],
                   [-sin_o*(omega**2),-cos_o*(omega**2),0],
                   [0,0,0]]) 
    
    zp[0:3] = z[3:6] 
    zp[3:6]= (-G*m_t/r**3)*z[0:3]  - RT@(R2@z[0:3] + 2*R1@z[3:6])
    return zp

# POSICION REAL SATELITE
t, x, y, z, vx, vy, vz = leer_eof('S1B_OPER_AUX_POEORB_OPOD_20200812T110757_V20200722T225942_20200724T005942.EOF') 

z0 = np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]])
zf = np.array([x[-1],y[-1],z[-1],vx[-1],vy[-1],vz[-1]])

# SOLUCIONES
sol = odeint(zpunto,z0,t)
x_sol = sol[:,0]
y_sol = sol[:,1]
z_sol = sol[:,2]

lista = []
norma_real = []
norma_sol = []
for i in range(9361):
    deriva = (np.sqrt((x[i]-x_sol[i])**2+(y[i]-y_sol[i])**2+(z[i]-z_sol[i])**2))/1000
    deriva1 = (np.sqrt((x[i])**2+(y[i])**2+(z[i])**2))/1000
    deriva2 = (np.sqrt((x_sol[i])**2+(y_sol[i])**2+(z_sol[i])**2))/1000
    lista.append(deriva)
    norma_real.append(deriva1)
    norma_sol.append(deriva2)
    
# GRAFICOS
x1=[0,18000,36000,54000,72000,90000]
xticks=["0","5","10","15","20","25"]
y1=[-5000000 , 0 , 5000000]
yticks=["-5000" , "0" , "5000"]

plt.figure()


# Grafico 1
plt.subplot(3,1,1)
plt.plot(t,x)
plt.plot(t,x_sol)
plt.xticks(x1,xticks)
plt.yticks(y1,yticks)
plt.ylabel("X [Km]")

plt.subplot(3,1,2)
plt.plot(t,y)
plt.plot(t,y_sol)
plt.xticks(x1,xticks)
plt.yticks(y1,yticks)
plt.ylabel("Y [Km]")

plt.subplot(3,1,3)
plt.plot(t,z)
plt.plot(t,z_sol)
plt.xticks(x1,xticks)
plt.yticks(y1,yticks)
plt.ylabel("Z [Km]")

plt.title("Posición en los ejes de coordenadas")
plt.xlabel("Tiempo, t [horas]")

plt.tight_layout()
plt.savefig('grafico_xyz.png')
plt.show()


# Grafico 2
plt.plot(t,lista)

plt.title("Distancia entre posicion real y predicha, δ = 597.9 [km]")
plt.xlabel('Tiempo, t [h]')
plt.ylabel("Deriva [Km]")
plt.xticks(x1,xticks)
plt.savefig('grafico_deriva.png')
plt.show()

