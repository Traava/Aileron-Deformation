import numpy as np
from variables import *

#------Loading the aerodynamic data file
aero_data = np.genfromtxt("aerodynamicloadf100.dat",delimiter = ",")

print(raw_data.shape)
#------Creating helper angles to find mesh coordinates
theta_x = np.zeros(Nx+1)
for j in range(Nx+1):
    theta_x[j] = (j)*np.pi/Nx

theta_z = np.zeros(Nz+1)
for i in range(Nz+1):
    theta_z[i] = (i)*np.pi/Nz

#-----Finding mesh coordinates
x_coor = np.zeros(Nx)
for i in range(Nx):
    x_coor[i] = (la*(1-np.cos(theta_x[i]))/2+la*(1-np.cos(theta_x[i+1]))/2)/2

z_coor = np.zeros(Nz)
for j in range(Nz):
    z_coor[j] = (Ca*(1-np.cos(theta_z[j]))/2+Ca*(1-np.cos(theta_z[j+1]))/2)/2

#print(x_coor[0],x_coor[-1], len(x_coor))
#------Lift force contribution for spanwise locations

q_tilde = np.zeros(Nx) #lift distribution value to be found for each x location

for i in range(Nx):
    chord_array = aero_data[:,i]
    chord_contrib = 0
    for j in range(len(chord_array)-1):
        chord_contrib += 






