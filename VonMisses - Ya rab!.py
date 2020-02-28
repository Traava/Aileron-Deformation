import numpy as np
from Shear_stress import *
from Bending_stress import *
from variables import *
import matplotlib.pyplot as plt


n    = 350 
n_sp = 100
m    = 60

x = np.linspace(0, la, m)

#ShearStress = tot_shear = np.hstack((tau1,tau3,tau4,tau6,tau5,tau2))
#DirectStress = Stress
max_vm = np.array([0,0])
for i in range(m):
    Sy = Vy(x)[i]
    Sz = Vz(x)[i]
    Ts  = T(x)[i]
    
    qb1 = ShearDistribution(Ts, Sy, Sz)[0]
    qb2 = ShearDistribution(Ts, Sy, Sz)[1]
    qb3 = ShearDistribution(Ts, Sy, Sz)[2]
    qb4 = ShearDistribution(Ts, Sy, Sz)[3]
    qb5 = ShearDistribution(Ts, Sy, Sz)[4]
    qb6 = ShearDistribution(Ts, Sy, Sz)[5]
    
    #qb2, qb3, qb4, qb5, qb6 = ShearDistribution(Ts, Sy, Sz)
    tau1 = qb1/t_sk
    tau2 = qb2/t_sp
    tau3 = qb3/t_sk
    tau4 = qb3/t_sk
    tau5 = qb5/t_sp
    tau6 = qb6/t_sk
    
    DirectStress = np.hstack((Bending(n, n_sp, m)[0][i],Bending(n, n_sp, m)[1][i]))
    ShearStress  = np.hstack((tau1,tau3,tau4,tau6,tau2,tau5))
    
    sigma_vm = np.sqrt(3*(ShearStress**2) + DirectStress**2)
    
    
    max_vm_x = np.array([i, np.max(sigma_vm)])
    max_vm   = np.vstack((max_vm, max_vm_x))
    


Max_location = np.where(max_vm == np.max(np.amax(max_vm, axis = 0)))
Max_loc = np.asscalar(Max_location[0])        
    
tau1 = ShearDistribution(T(x)[Max_loc], Vy(x)[Max_loc], Vz(x)[Max_loc])[0]/t_sk
tau2 = ShearDistribution(T(x)[Max_loc], Vy(x)[Max_loc], Vz(x)[Max_loc])[1]/t_sp
tau3 = ShearDistribution(T(x)[Max_loc], Vy(x)[Max_loc], Vz(x)[Max_loc])[2]/t_sk
tau4 = ShearDistribution(T(x)[Max_loc], Vy(x)[Max_loc], Vz(x)[Max_loc])[3]/t_sk
tau5 = ShearDistribution(T(x)[Max_loc], Vy(x)[Max_loc], Vz(x)[Max_loc])[4]/t_sp
tau6 = ShearDistribution(T(x)[Max_loc], Vy(x)[Max_loc], Vz(x)[Max_loc])[5]/t_sk

DirectStress = np.hstack((Bending(n, n_sp, m)[0][Max_loc],Bending(n, n_sp, m)[1][Max_loc]))
ShearStress  = np.hstack((tau1,tau3,tau4,tau6,tau2,tau5))
#ShearStress = np.zeros(len(Direct_stress))

sigma_vm = np.sqrt(3*(ShearStress**2) + DirectStress**2)

fig1 = plt.figure()
ax1 = plt.axes()
ax1.invert_yaxis()
img1 = plt.scatter(Bending(n, n_sp, m)[2], Bending(n, n_sp, m)[3], c = sigma_vm, cmap ='RdYlGn_r')
#img2 = plt.scatter(z_coord_sp, y_coord_sp, c = Stress_sp[20], cmap ='RdYlGn_r') 


cbar1 = fig1.colorbar(img1)
cbar1.set_label(r'[$\frac{N}{m^2}$]')
ax1.set_title("Von Mises Stress Distribution")
ax1.set_xlabel('z [m]')
ax1.set_ylabel('y [m]')











