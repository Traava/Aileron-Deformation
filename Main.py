import numpy as np
from numpy import sin, cos
import scipy as sp
from variables import *

#test variables:

Izz = 2
G = 9 
J = 9

E = 10e7
theta  = 30
zh = 0.1
hh = zh - h/2


A = np.zeros((9,9))


def Mac(x):
    return (max(x,0))

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
F1_const = 1/E*Izz*Q4 + hh/(G*J)*tao2

#constant

################################################################

#########################Form 2################################## w(x)
def F2_Rz123(x,xn):
    return 1/(6*E*Iyy)*Mac(x-xn)**3

def F2_RaIaII(x,xn):
    return 1/(6*E*Iyy)*cos(theta)*Mac(x-xn)**3


################################################################
#########################Form 3A#################################### ( v(x)+zh*phi(x) )* sin(theta)

def F3A_Ry123(x,xn):
    return ( 1/(6*E*Izz)*Mac(x-xn)**3 + zh*hh/(G*J)*Mac(x-xn) )*sin(theta)

def F3A_RaIaII(x, xn):
    return ( 1/(6*E*Izz)* sin(theta) * Mac(x - xn)**3 +zh/(G*J)*Mac(x-xn)*(zh*sin(theta)-h/2*cos(theta)) )*sin(theta)

F3A_const = ( 1/E*Izz*Q4 + zh/(G*J)*tao2 )*sin(theta)

##############################################################################

#########################Form 3B##################################### ( w(x) -h/2 *phi(x) )

def F3B_Ry123(x,xn):
    return  -(  h*hh/(2*G*J)*Mac(x-xn) )*cos(theta)
def F3B_Rz123(x,xn):
    return ( 1/(6*E*Izz)*Mac(x-xn)**3 )*cos(theta)

def F3B_RaIaII(x, xn):
    return ( 1/(6*E*Izz)* cos(theta) * Mac(x - xn)**3 -h/(2*G*J)*Mac(x-xn)*(zh*sin(theta)-h/2*cos(theta)) )*cos(theta)

F3B_const = -h/(2*G*J)* tao2 *cos(theta)

##################################################################

#IMPLENTING THE MATRIX

 X = [Ry1,Ry2,Ry3,Rz1,Rz2,Rz3,RaI,C1,C2,C3,C4,C5]  

a1 = [1,1,1,0,0,0,sin(theta), 0,0,0,0,0]
a2 = [0,0,0,1,1,1,cos(theta), 0,0,0,0,0] 
a3 = [0,0,0,0,0,0,h*(sin(theta)-cos(theta))/2,0,0,0,0,0]
a4 = [0,0,0,0, x1-x2, x1-x3, (x1-xai)*cos(theta), 0 , 0,0,0,0]
a5 = [xai -x1, x2-xai, x3-xai, 0,0,0,0,0,0,0,0,0]
a6 = [F1_Ry123(x1,x1), F1_Ry123(x1,x2), F1_Ry123(x1, x3), 0,0,0,F1_RaIaII(x1, xaI), x1, 1,0,0,hh]
a7 = [0,0,0, F2_Rz123(x1,x1), F2_Rz123(x1,x2), F2_Rz123(x1,x3), F2_RaIaII(x1, xaI), 0, 0, x1, 1, 0]
a8 = [F1_Ry123(x2,x1), F1_Ry123(x2,x2), F1_Ry123(x2, x3), 0,0,0,F1_RaIaII(x2, xaI), x2, 1,0,0,hh]
a9 = [0,0,0, F2_Rz123(x2,x1), F2_Rz123(x2,x2), F2_Rz123(x2,x3), F2_RaIaII(x2, xaI), 0, 0, x2, 1, 0]
a10 = [F1_Ry123(x3,x1), F1_Ry123(x3,x2), F1_Ry123(x3, x3), 0,0,0,F1_RaIaII(x3, xaI), x3, 1,0,0,hh]
a11 = [0,0,0, F2_Rz123(x3,x1), F2_Rz123(x3,x2), F2_Rz123(x3,x3), F2_RaIaII(x3, xaI), 0, 0, x3, 1, 0]
a12 = [F3A_Ry123(xaI, x1)+F3B_Ry123(xaI, x1), F3A_Ry123(xaI,x2)+F3B_Ry123(xaI,x2), F3A_Ry123(xaI, x3)+F3B_Ry123(xaI, x3), 
       F3B_Rz123(xaI, x1), F3B_Rz123(xaI, x2), F3B_Rz123(xaI, x3), F3A_RaIaII(xaI,xaI)+F3B_RaIaII(xaI,xaI), xaI, 1, 0,0,hh]

 
b1 = [-P*sin(theta) + Q1]
b2 = [-P*sin(theta)]
b3 = [P*h/2*(cos(theta)-sin(theta)) + Q1 * ]##
b4 = []
b5 = []
b6 = [-d1*cos(theta) - F1_RaIaII(x1,xaII) - 1/(E*Izz)*Q4 - hh/(G*J)*tao2]
b7 = []
b8 = [- F1_RaIaII(x2,xaII) - 1/(E*Izz)*Q4 - hh/(G*J)*tao2]
b9 =[]
b10 = []
b11 = []
b12 = []
 
 
 
 
 
 
 
 