#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 05:07:48 2020

@author: Tomas
"""

import  matplotlib.pylab as plt
import numpy as np
from scipy import pi , sin , cos , dot , sqrt , zeros
from leereof import leer_eof
from scipy.integrate import odeint
from time import perf_counter

# COMIENZA EL CODIGO
t1 = perf_counter()

# PARAMETROS
G = 6.67*10**-11 
m_t = 5.98*10**24 
omega = 2.*pi/86400 

# FUNCIONES
def zpunto(z,t):
    #parámetros
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
    
    J2 = 1.75553*(10**10)*1000**5
    J3 = -2.61913*(10**11)*1000**6
    x = z[0]
    y = z[1]
    zz = z[2]
    
    #defino mis valores de F
    der_x_J2 = J2*(x/r**7)*(6*zz**2-3/2*(x**2+y**2))
    der_y_J2 = J2*(y/r**7)*(6*zz**2-3/2*(x**2+y**2))
    der_z_J2 = J2*(zz/r**7)*(3*zz**2-9/2*(x**2+y**2))
    FJ2 = np.array([der_x_J2,der_y_J2,der_z_J2])
    
    der_x_J3 = J3*x*zz/r**9*(10*zz**2-15/2*(x**2+y**2))
    der_y_J3 = J3*y*zz/r**9*(10*zz**2-15/2*(x**2+y**2))
    der_z_J3 = J3/r**9*(4*zz**2*(zz**2-3*(x**2+y**2))+3/2*(x**2+y**2)**2)
    FJ3 = np.array([der_x_J3,der_y_J3,der_z_J3])
    
    zp[0:3] = z[3:6]
    z2 = ((-G*m_t)/(r**3))*z[0:3] - RT@(R2@z[0:3] + 2*R1@z[3:6]) - FJ2 + FJ3
    zp[3:6] = z2
    return zp

def eulerint(zp,z0,t,Nsub):
    Nt = len(t)
    Ndim = len(np.array(z0))
    z = np.zeros((Nt,Ndim))
    z[0,:] = z0[:]
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsub
        z_temp = z[i-1,:].copy()
        for k in range(Nsub):
            z_temp += dt * zpunto(z_temp,t_anterior + k*dt)
        z[i,:] = z_temp
    return z

# POSICION REAL SATELITE
t, x, y, z, vx, vy, vz = leer_eof('S1B_OPER_AUX_POEORB_OPOD_20200812T110757_V20200722T225942_20200724T005942.EOF') 
z0 = np.array([x[0],y[0],z[0],vx[0],vy[0],vz[0]]) 

# SOLUCIONES
sol = odeint(zpunto,z0,t)


x_sol = sol[:,0]
y_sol = sol[:,1]
z_sol = sol[:,2]

x_ode = sol[:,0]
y_ode = sol[:,1]
z_ode = sol[:,2]

zp = sol[:,:] 
t3 = perf_counter()
sol2 = eulerint(zp,z0,t,1)
t4 = perf_counter()

x_euler = sol2[:,0]
y_euler = sol2[:,1]
z_euler = sol2[:,2]


# DERIVA FINAL
delta = np.sqrt((x-x_ode)**2+(y-y_ode)**2+(z-z_ode)**2)
delta2 = np.sqrt((x-x_euler)**2+(y-y_euler)**2+(z-z_euler)**2)

delta_final = (np.sqrt((x_ode[-1]-x_euler[-1])**2+(y_ode[-1]-y_euler[-1])**2+(z_ode[-1]-z_euler[-1])**2))/1000


# LISTA DISTANCIAS
lista = []
for i in range(9361):
    deriva = (np.sqrt((x[i]-x_sol[i])**2+(y[i]-y_sol[i])**2+(z[i]-z_sol[i])**2))/1000
    lista.append(deriva)

print (lista[-1])

# GRAFICOS

# Grafico 1
x1=[0,18000,36000,54000,72000,90000]
xticks=["0","5","10","15","20","25"]
y1=[-5000000 , 0 , 5000000]
yticks=["-5000" , "0" , "5000"]

plt.figure()

plt.subplot(3,1,1)
plt.title("Posición en los ejes de coordenadas con J2 y J3")
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

plt.xlabel("Tiempo, t [horas]")

plt.tight_layout()
plt.savefig('grafico_xyz_J2yJ3.png')
plt.show()


# Grafico 2
plt.plot(t,lista)

plt.title("Distancia entre posicion real y predicha con J2 y J3, δ = 317.59 [km]")
plt.xlabel('Tiempo, t [h]')
plt.ylabel("Deriva [Km]")
plt.xticks(x1,xticks)
plt.savefig('grafico_deriva_J2yJ3.png')
plt.show()


# Grafico 3
plt.plot(t , delta , "r" , label = "Odeint")
plt.plot(t , delta2 , "b" , label = "Eulerint")

x1=[ 0 , 18000 , 36000 , 54000 , 72000 , 90000]
xticks=["0" , "5" , "10" , "15" , "20" , "25"]
plt.xticks(x1,xticks)

plt.title("Gráfico de deriva con J2 y J3, Nsub = 1")
plt.ylabel('Delta [km]')
plt.xlabel('Tiempo, t [horas]')
plt.tight_layout()
plt.grid(True)
plt.legend()

plt.savefig('ode_vs_eul_deriva_nsub=1_J2yJ3')
plt.show()

t2 = perf_counter()
dif_ode = t2 - t1
print (dif_ode)
