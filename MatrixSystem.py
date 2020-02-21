import numpy as np
from numpy import sin, cos
import scipy as sp
from variables import *
import Aero_load as al
import MoI



#test variables:

Izz = MoI.Izz_Aileron
Iyy = MoI.Iyy_Aileron
J = Iyy+Izz
theta = theta_max
zh = zhat
hh = zh- h/2



A = np.zeros((9,9))

#Macaulay step function
def Mac(x):
    return (max(x,0))
#INTEGRALS
def Q1(x):
    return al.integration(1, x, al.q_tilde, al.x_coor)

def Q4(x):
    return al.integration(4, x, al.q_tilde, al.x_coor)
def tao2(x):
    return al.integration(4, x, al.qT_tilde, al.x_coor)
def aero_Mx(x):
    return al.integration(1,x,(al.CoPs-np.ones(len(al.CoPs))*h/2)*al.q_tilde, al.x_coor)
def aero_Mz(x):
    return al.integration(1,x, al.q_tilde*(x-xaI), al.x_coor)

#########################Form 1####################### Form 1: v(x) +hh*phi(x)  

#coeff. for Ry1, Ry2, Ry3:
def F1_Ry123(x, xn):
    return 1/(6*E*Izz)*Mac(x-xn)**3 + hh/(G*J)*Mac(x-xn)

#coeff. for RaI, RaII:
def F1_RaIaII(x, xn):
    return 1/(6*E*Izz)* sin(theta) * Mac(x - xn)**3 +hh/(G*J)*Mac(x-xn)*(zh*sin(theta)-h/2*cos(theta))

#coeff for C1
def F1_C1(x):
    return x

#coeff for C2
F1_C2 = 1

#constant values 
#Q4 is four-time-integral of q(x)
#tao2 is torque of aerodynamics

def F1_const(x):
    return 1/E*Izz*Q4(x)  + hh/(G*J)* tao2(x)

#constant

################################################################

#########################Form 2################################## w(x)
def F2_Rz123(x,xn):
    return 1/(6*E*Iyy)*Mac(x-xn)**3

def F2_RaIaII(x,xn):
    return 1/(6*E*Iyy)*cos(theta)*Mac(x-xn)**3

def F2_const(x):
    return 1/(6*E*Iyy)*P*cos(theta)*Mac(x-xaII)**3


################################################################
#########################Form 3A#################################### ( v(x)+zh*phi(x) )* sin(theta)

def F3A_Ry123(x,xn):
    return ( 1/(6*E*Izz)*Mac(x-xn)**3 + zh*hh/(G*J)*Mac(x-xn) )*sin(theta)

def F3A_RaIaII(x, xn):
    return ( 1/(6*E*Izz)* sin(theta) * Mac(x - xn)**3 +zh/(G*J)*Mac(x-xn)*(zh*sin(theta)-h/2*cos(theta)) )*sin(theta)

def F3A_const(x):
    return ( 1/E*Izz*Q4(x) + zh/(G*J)*tao2(x) )*sin(theta)

##############################################################################

#########################Form 3B##################################### ( w(x) -h/2 *phi(x) )

def F3B_Ry123(x,xn):
    return  -(  h*hh/(2*G*J)*Mac(x-xn) )*cos(theta)
def F3B_Rz123(x,xn):
    return ( 1/(6*E*Izz)*Mac(x-xn)**3 )*cos(theta)

def F3B_RaIaII(x, xn):
    return ( 1/(6*E*Izz)* cos(theta) * Mac(x - xn)**3 -h/(2*G*J)*Mac(x-xn)*(zh*sin(theta)-h/2*cos(theta)) )*cos(theta)

def F3B_const(x):
    return -h/(2*G*J)* tao2(x) *cos(theta)

##################################################################

#IMPLENTING THE MATRIX

Xdef = ['Ry1','Ry2','Ry3','Rz1','Rz2','Rz3','RaI','C1','C2','C3','C4','C5']  

a1 = [1,1,1,0,0,0,sin(theta), 0,0,0,0,0]
a2 = [0,0,0,1,1,1,cos(theta), 0,0,0,0,0] 
a3 = [0,0,0,0,0,0,h*(sin(theta)-cos(theta))/2,0,0,0,0,0]
a4 = [0,0,0,0, x1-x2, x1-x3, (x1-xaI)*cos(theta), 0 , 0,0,0,0]
a5 = [xaI -x1, x2-xaI, x3-xaI, 0,0,0,0,0,0,0,0,0]
a6 = [F1_Ry123(x1,x1), F1_Ry123(x1,x2), F1_Ry123(x1, x3), 0,0,0,F1_RaIaII(x1, xaI), x1, 1,0,0,hh]
a7 = [0,0,0, F2_Rz123(x1,x1), F2_Rz123(x1,x2), F2_Rz123(x1,x3), F2_RaIaII(x1, xaI), 0, 0, x1, 1, 0]
a8 = [F1_Ry123(x2,x1), F1_Ry123(x2,x2), F1_Ry123(x2, x3), 0,0,0,F1_RaIaII(x2, xaI), x2, 1,0,0,hh]
a9 = [0,0,0, F2_Rz123(x2,x1), F2_Rz123(x2,x2), F2_Rz123(x2,x3), F2_RaIaII(x2, xaI), 0, 0, x2, 1, 0]
a10 = [F1_Ry123(x3,x1), F1_Ry123(x3,x2), F1_Ry123(x3, x3), 0,0,0,F1_RaIaII(x3, xaI), x3, 1,0,0,hh]
a11 = [0,0,0, F2_Rz123(x3,x1), F2_Rz123(x3,x2), F2_Rz123(x3,x3), F2_RaIaII(x3, xaI), 0, 0, x3, 1, 0]
a12 = [F3A_Ry123(xaI, x1)+F3B_Ry123(xaI, x1), F3A_Ry123(xaI,x2)+F3B_Ry123(xaI,x2), F3A_Ry123(xaI, x3)+F3B_Ry123(xaI, x3), 
       F3B_Rz123(xaI, x1), F3B_Rz123(xaI, x2), F3B_Rz123(xaI, x3), F3A_RaIaII(xaI,xaI)+F3B_RaIaII(xaI,xaI), xaI, 1, 0,0,hh]

 
b1 = [-P*sin(theta) + Q1(la-0.01)]
b2 = [-P*sin(theta)]
b3 = [P*h/2*(cos(theta)-sin(theta)) + aero_Mx(la-0.01) ]
b4 = [P*cos(theta)*(xaII-x1)]
b5 = [-P*sin(theta)*(xaII-xaI) - aero_Mz(la-0.01)]
b6 = [-d1*cos(theta) - F1_const(x1)]
b7 = [d1*sin(theta) - F2_const(x1)]
b8 = [- F1_const(x2)]
b9 = [- F2_const(x2)]
b10 = [-d3*cos(theta) - F1_const(x3)]
b11 = [d3*sin(theta) - F2_const(x3)]
b12 = [-F3A_const(xaI) - F3B_const(xaI)]

alst = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12]

blst = [b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12]

A = np.zeros((12,12))
B = np.zeros((12,1))

for i in range(len(alst)):
    assert len(alst[i]) == 12
    assert len(blst[i]) == 1

    A[i] = np.array(alst[i])
    B[i] = np.array(blst[i])

print()
print('Start Main')
print()
X = np.linalg.solve(A,B)

for i in range(len(X)):
    print(Xdef[i], '\t=\t ', X[i][0])
    
Ry1 = X[0]
Ry2 = X[1]
Ry3 = X[2]
Rz1 = X[3] 
Rz2 = X[4]
Rz3 = X[5]
RaI = X[6]
C1 = X[7]
C2 = X[8]
C3 = X[9]
C4 = X[10]
C5  = X[11]
    
sumFy = Ry1 + Ry2 + Ry3 + RaI*sin(theta) + P*sin(theta) - Q1(la-0.01)
 
 