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


# ------- Cubic spline interpolation

def interpol(f,x,M0,MN):
    dim = len(x)
    system_mat = np.zeros((dim,dim))            #System looks like: system_mat*M = output_vec
    system_mat[0,0] = 1                         #to say that M0 = M0
    system_mat[-1,-1] = 1                       #to say that MN = MN
    output_vec = np.zeros(dim)
    output_vec[0] = M0                          #to say that M0 = M0
    output_vec[-1]= MN                          #to say that MN = MN
    h=np.zeros(dim-1)
    h[0] = x[1]-x[0]
    for i in range(1,dim-1):                    #Check Numerical Analysis Reader page 43
        h_prev = x[i]  - x[i-1]
        h_i    = x[i+1]- x[i]
        h[i]   = h_i
        system_mat[i,i-1] = h_prev/6
        system_mat[i,i]   = (h_prev + h_i)/3
        system_mat[i,i+1] = h_i/6
        output_vec[i] = (f[i+1] - f[i])/h_i - (f[i]-f[i-1])/h_prev

    M = np.linalg.solve(system_mat, output_vec)
    coef_mat = np.zeros((dim-1,4))
    for i in range(dim-1):
        print(h[i])
        ai = (M[i+1]-M[i])/(6*h[i])
        bi = M[i]/2
        ci = (f[i+1]-f[i])/h[i]-h[i]*M[i]/3 - h[i]*M[i+1]/6
        di = f[i]
        Di = di - ci*x[i] + bi*(x[i])**2 -ai*(x[i])**3
        Ci = ci - 2*x[i]*bi + 3*ai*(x[i])**2
        Bi = bi - 3*ai*x[i]
        Ai = ai
        coef_mat[i] = np.array([Ai,Bi,Ci,Di])
    return coef_mat



coef_mat = interpol(q_tilde,x_coor,0,0)

def polynomial(x,A,B,C,D):
    return A*x**3 + B*x**2 + C*x+D

points = []
values = []
n = 10
for i in range(len(x_coor)-1):
    h_i = x_coor[i + 1] - x_coor[i]
    A ,B, C, D = tuple(coef_mat[i])
    for j in range(n):
        points.append(x_coor[i]+h_i*j/n)
        values.append(polynomial(x_coor[i]+h_i*j/n,A,B,C,D))

def int_pol(coefs_in):
    size = len(coefs_in)
    syst_mat = np.zeros((size+1,size))
    for i in range(1,size+1):
        syst_mat[i,i-1] = 1/i
    coefs_out = np.dot(syst_mat,coefs_in)
    return coefs_out

print(int_pol(int_pol(np.array([1,1,1,1]).reshape(-1,1))))




