
import  matplotlib.pylab as plt
import numpy as np
from scipy import pi , sin , cos , dot , sqrt , zeros
from leereof import leer_eof
from scipy.integrate import odeint
from time import perf_counter

# PARAMETROS
G = 6.67 * (10**-11) 
m_t = 5.98 * (10**24) 
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
     
    zp[0:3] = z[3:6] 
    zp[3:6]= (-G*m_t/r**3)*z[0:3]  - RT@(R2@z[0:3] + 2*R1@z[3:6])
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
t1 = perf_counter()
sol = odeint(zpunto,z0,t)
t2 = perf_counter()
x_ode = sol[:,0]
y_ode = sol[:,1]
z_ode = sol[:,2]

dif_ode = t2 - t1
print (dif_ode)

zp = sol[:,:] 
t3 = perf_counter()
sol2 = eulerint(zp,z0,t,1700)
t4 = perf_counter()

dif_eul = t4 - t3
print (dif_eul)

x_euler = sol2[:,0]
y_euler = sol2[:,1]
z_euler = sol2[:,2]

# DERIVA FINAL
delta = np.sqrt((x-x_ode)**2+(y-y_ode)**2+(z-z_ode)**2)
delta2 = np.sqrt((x-x_euler)**2+(y-y_euler)**2+(z-z_euler)**2)

delta_final = (np.sqrt((x_ode[-1]-x_euler[-1])**2+(y_ode[-1]-y_euler[-1])**2+(z_ode[-1]-z_euler[-1])**2))/1000


# GRAFICO
plt.plot(t , delta , "r" , label = "Odeint")
plt.plot(t , delta2 , "b" , label = "Eulerint")

x1=[ 0 , 18000 , 36000 , 54000 , 72000 , 90000]
xticks=["0" , "5" , "10" , "15" , "20" , "25"]
y1=[ 0 , 5000000 , 10000000 , 15000000 , 20000000 ]
yticks=[ "0" , "5000" , "10000" , "15000" , "20000"]
plt.xticks(x1,xticks)
plt.yticks(y1,yticks)

plt.title("Gráfico de deriva , Nsub = 1700")
plt.ylabel('Delta [km]')
plt.xlabel('Tiempo, t [horas]')
plt.tight_layout()
plt.grid(True)
plt.legend()

plt.savefig('ode_vs_eul_deriva_nsub=1700')
plt.show()
