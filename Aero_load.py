import numpy as np
from variables import *
import matplotlib.pyplot as plt

#------Loading the aerodynamic data file
aero_data = np.genfromtxt("aerodynamicloadf100.dat",delimiter = ",")

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


#----- Getting spanwise lift distribution
q_tilde  = np.zeros(Nx)                                                     #lift distribution value to be found for each x location
qT_tilde = np.zeros(Nx)                                                     #torque of aero load for each x loc taken around leading edge
CoPs     = np.zeros(Nx)                                                     #center of pressure for each x loc

for i in range(Nx):
    chord_array    = aero_data[:,i]                                         #picking out lift data for given spanwise location
    lift_contrib   = 0                                                      #the integral of lift over that spanwise location
    torque_contrib = 0                                                      #integral of torque over that spanwise location taken around leading edge
    for k in range(len(chord_array)-1):
        segm_len        = z_coor[k+1]-z_coor[k]                                    #distance between two adjacent chord locations
        dA              = segm_len*(chord_array[k]+chord_array[k+1])/2             #trapezoidal rule
        lift_contrib   += dA                                                       #lift sums areas
        torque_contrib += (z_coor[k]+segm_len/2)*dA                                #torque sums weighted areas
    q_tilde[i]  = lift_contrib
    qT_tilde[i] = torque_contrib
    CoPs[i]     = torque_contrib/lift_contrib                               #Distance is moment/force







