import numpy as np
from variables import *
import matplotlib.pyplot as plt

"UNIT 2.1"
#------Loading the aerodynamic data file
aero_data = 1000*np.genfromtxt("aerodynamicloadf100.dat",delimiter = ",")  #1000* to change to N/m**2

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

"UNIT 2.2"
#----- Getting spanwise lift distribution
q_tilde  = np.zeros(Nx)          #lift distribution value to be found for each x location
qT_tilde = np.zeros(Nx)          #torque of aero load for each x loc taken around shear center
CoPs     = np.zeros(Nx)          #center of pressure for each x loc

for i in range(Nx):
    chord_array    = aero_data[:,i]     #picking out lift data for given spanwise location
    lift_contrib   = 0                  #the integral of lift over that spanwise location
    torque_contrib = 0                  #integral of torque over that spanwise location taken around shear center
    torque_le      = 0

    for k in range(len(chord_array)-1):
        segm_len        = z_coor[k+1]-z_coor[k]                         #distance between two adjacent chord locations
        dA              = segm_len*(chord_array[k]+chord_array[k+1])/2  #trapezoidal rule
        lift_contrib   += dA                                            #lift sums areas
        torque_le      += (z_coor[k]+segm_len/2)*dA                     #torque sums weighted areas
        torque_contrib += (zhat - z_coor[k]-segm_len/2)*dA              #torque sums weighted areas

    q_tilde[i]  = lift_contrib
    qT_tilde[i] = torque_contrib
    CoPs[i]     = torque_le/lift_contrib                               #Distance is moment/force

"UNIT 2.3"
# ------- Cubic spline interpolation
def interpol(f,coor,M0,MN):

    dim = len(coor)
    system_mat = np.zeros((dim,dim))            #System looks like: system_mat*M = output_vec
    system_mat[0,0] = 1                         #to say that M0 = M0
    system_mat[-1,-1] = 1                       #to say that MN = MN
    output_vec = np.zeros(dim)
    output_vec[0] = M0                          #to say that M0 = M0
    output_vec[-1]= MN                          #to say that MN = MN
    h=np.zeros(dim-1)
    h[0] = coor[1]-coor[0]

    for i in range(1,dim-1):                    #Check Numerical Analysis Reader page 43
        h_prev            = coor[i]  - coor[i-1]
        h_i               = coor[i+1]- coor[i]
        h[i]              = h_i
        system_mat[i,i-1] = h_prev/6
        system_mat[i,i]   = (h_prev + h_i)/3
        system_mat[i,i+1] = h_i/6
        output_vec[i]     = (f[i+1] - f[i])/h_i - (f[i]-f[i-1])/h_prev

    M = np.linalg.solve(system_mat, output_vec)#to find the M coefficients
    coef_mat = np.zeros((dim-1,4))             #In between every two neighbouring x positions we have a spline (dim-1 of them), each spline has 4 coefs Ax3 + Bx2 + Cx +D

    for i in range(dim-1):
        # numerical analysis again
        ai          = (M[i+1]-M[i])/(6*h[i])
        bi          = M[i]/2
        ci          = (f[i+1]-f[i])/h[i]-h[i]*M[i]/3 - h[i]*M[i+1]/6
        di          = f[i]
        # Reformating to get the desired equation in form Ax3 + Bx2 + Cx +D
        Di          = di - ci*coor[i] + bi*(coor[i])**2 -ai*(coor[i])**3
        Ci          = ci - 2*coor[i]*bi + 3*ai*(coor[i])**2
        Bi          = bi - 3*ai*coor[i]
        Ai          = ai
        coef_mat[i] = np.array([Di,Ci,Bi,Ai])

    return coef_mat

"UNIT 2.4"
#Evaluating a polynomial at one point
def polynomial(x,coefs):

    result = 0
    for k in range(len(coefs)):
        result+= coefs[k]*x**(k)         #the kth power corresponds to the kth coefficient

    return result

"UNIT 2.5"
#------Function that interpolates polynpomials analytically and gives definite integral value
def int_pol(coefs_in,a,b):  #coefs are A,B,C,D etc. coefficients of a polynom. given as an array, a and b are start and end of definite integral resp.

    size = len(coefs_in)
    syst_mat = np.zeros((size+1,size))

    for i in range(1,size+1):
        syst_mat[i,i-1] = 1/i                                   #this is explained more in depth in the report. Esentially creating an integration in the form of linear transformation
    coefs_out = np.dot(syst_mat,coefs_in)
    result = polynomial(b,coefs_out)-polynomial(a,coefs_out)    #finally evaluating F(b) - F(a)

    return result

"UNIT 2.6"
#-------Function that iterates over the span and finds function values for the antiderivative
def int_span(coef_mat,coor):

    dim = len(coor)
    contributions = np.zeros(dim-1)

    for i in range(dim-1):
        contributions[i] = int_pol(coef_mat[i],coor[i],coor[i+1])  #each spanwise segment has an area asociated with it

    integrated_vals = np.zeros(dim)

    for j in range(1,dim):
        integrated = 0
        for k in range(j):                                         #need to sum all contributions up until the point we are looking at
            integrated += contributions[k-1]
        integrated_vals[j] = integrated

    return integrated_vals

"UNIT 2.7"
#-------Function which combines the two above, first integrating and then interpolating again
def int_ana(coef_mat,coor):

    new_mat = interpol(int_span(coef_mat,coor),coor,0,0)

    return new_mat

"UNIT 2.8"
#-------Finally, a function to which we give point x and a given number of integrations and it evaluates this for us
def integration(n,x,func,coor):

    func_coefs = interpol(func,coor,0,0)                    #first interpolant

    for i in range(n):                                      #integrate n times
        func_coefs = int_ana(func_coefs,coor)

    if x<coor[0]:
        return 0
    else:
        for j in range(1,len(coor)):
            if coor[j] >= x:                                     #this to find the actual value corresponding to our x
                output = polynomial(x,func_coefs[j-1])
                break

    return output




#VERIFICATION BELOW
#print(CoPs/Ca)
la = x_coor[-1]
print(integration(1,la,q_tilde,x_coor))
print(x_coor)
print(la)


"UNIT 2.1 Verification"

x_cr = []
z_cr = []

#Creating the mesh coordinates (xi,zj)
for xi in x_coor:
    for zj in z_coor:
        x_cr.append(xi)
        z_cr.append(zj)

#Plotting
plt.scatter(x_cr,z_cr,s = 3)
plt.ylim(ymin=0,ymax = Ca)
plt.xlim(xmin = 0, xmax = la)
plt.ylabel("z [m]").set_size(15)
plt.xlabel("x [m]").set_size(15)
plt.show()

#Value verification

print("Smallest value of x is: ",x_coor[0])
print("Smallest value of z is: ",z_coor[0])
print("Largest value of x is: ",x_coor[-1])
print("Largest value of z is: ",z_coor[-1])

#----- This is basically for visual verification of the preciseness of the spline interpolation
# coef_mat = interpol(q_tilde,x_coor,0,0)
# points = []
# values = []
# n = 10
# for i in range(len(x_coor)-1):
#     h_i = x_coor[i + 1] - x_coor[i]
#     coefs = coef_mat[i]
#     for j in range(n):
#         points.append(x_coor[i]+h_i*j/n)
#         values.append(polynomial(x_coor[i]+h_i*j/n,coefs))
#
# plt.plot(points,values)
# plt.show()

print(integration(1,la,q_tilde,x_coor))
print(x_coor)