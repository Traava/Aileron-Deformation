import numpy as np
from numpy import sin, cos
import scipy as sp
from variables import *
import Aero_load as al
import MoI
import seaborn as sns
sns.set()


#variables:

Izz = MoI.Izz_Aileron
Iyy = MoI.Iyy_Aileron
J =7.748548555816593e-06     #0.000175317141188
theta = theta_max
zh = 0.08553893540215983    #zhat
hh = zh- h/2

la = al.x_coor[-1]
mult = 1
A = np.zeros((9,9))

print(MoI.Izz_Aileron)
print(MoI.Iyy_Aileron)
"Unit 3.1"
#Macaulay step function
def Mac(x):
    if np.ndim(x) == 0:
        return (max(x,0))
    elif np.ndim(x) == 1:
        negatives = x < 0
        x[negatives] = 0
        return x
    else:
        return 'Failed Macaulay'
    
def MacTest():
    #array of test values
    testVal = np.arange(-8, 12, 0.97)
    
    #testing Macaulay of scalars
    for val in testVal:
        print('[',val,'] = ',Mac(val))
    print()
    #testing Macaulay of array
    print(testVal)
    print(Mac(testVal))
     
    
"Unit 3.2"    
#Integral of aerodynamic load distribution over x
def Q1(x):
    if np.ndim(x) == 0:
        return al.integration(1, x, mult*al.q_tilde, al.x_coor)
    elif np.ndim(x) == 1:
        integ = np.zeros(np.shape(x))
        for i in x:
            integ[i] = al.integration(1, i, mult*al.q_tilde, al.x_coor)
        return integ

#four-time integral of aerodynamic load distribution over x
def Q4(x):
    if np.ndim(x) == 0:
        return al.integration(4, x, mult*al.q_tilde, al.x_coor)
    elif np.ndim(x) == 1:
        integ = np.zeros(np.shape(x))
        for i in range(len(x)):
            integ[i] = al.integration(4, x[i], mult*al.q_tilde, al.x_coor)
        return integ
    
#Q_test = []
#for x in al.x_coor:
#    Q_test.append(Q1(x))
#    
#plt.plot(al.x_coor, Q_test, al.x_coor, al.q_tilde)
#plt.show()
    

#Two-time integral of torque distribution over x"
def tao2(x):
    if np.ndim(x) == 0:
        return al.integration(2, x, mult*al.qT_tilde, al.x_coor)
    elif np.ndim(x) == 1:
        integ = np.zeros(np.shape(x))
        for i in range(len(x)):
            integ[i] = al.integration(2, x[i], mult*al.qT_tilde, al.x_coor)
        return integ
    
"Unit 3.3"
def aero_Mx(x):
    return al.integration(1,x,-(al.CoPs-np.ones(len(al.CoPs))*h/2)*mult*al.q_tilde, al.x_coor)       # Marek edit: added minus
def aero_Mz(x):
    return al.integration(1,x, mult*al.q_tilde*(al.x_coor-np.ones(len(al.x_coor))*xaI), al.x_coor)   #Marek edit: the arm was wrong here

#########################Form 1####################### Form 1: v(x) +hh*phi(x)  
print("Hello this is me:", al.integration(1,la,mult*al.q_tilde*(al.x_coor-np.ones(len(al.x_coor))*xaI), al.x_coor))

"Unit 3.4"
#coeff. for Ry1, Ry2, Ry3:
def F1_Ry123(x, xn):
    return 1/(6*E*Izz)*Mac(x-xn)**3 + (hh)**2/(G*J)*Mac(x-xn) 

#coeff. for RaI, RaII:
def F1_RaIaII(x, xn):
    return 1/(6*E*Izz)* sin(theta) * Mac(x - xn)**3 +hh/(G*J)*Mac(x-xn)*(zh*sin(theta)-h/2*cos(theta))


#constant values 
def F1_const(x):
    return 1/(E*Izz)*Q4(x)  + hh/(G*J)* tao2(x)



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

def F3P_const(x,xn):
    return P*(((1/(E*Izz*6))*(sin(theta))**2 + (1/(6*Iyy*E))*(cos(theta))**2)*Mac(x-xn)**3 + ((zh*sin(theta)-h*cos(theta)/2)/(G*J))*(sin(theta)*zh-cos(theta)*h/2)*Mac(x-xn))

##################################################################

"Unit 3.5"
#IMPLENTING THE MATRIX
Xdef = ['Ry1','Ry2','Ry3','Rz1','Rz2','Rz3','RaI','C1','C2','C3','C4','C5']  

a1 = [1,1,1,0,0,0,sin(theta), 0,0,0,0,0]
a2 = [0,0,0,1,1,1,cos(theta), 0,0,0,0,0] 
a3 = [hh,hh,hh,0,0,0,zh*sin(theta)-h*cos(theta)/2,0,0,0,0,0]     #could try torque
a4 = [0,0,0,x1-la, x2-la, x3-la, (xaI-la)*cos(theta), 0 , 0,0,0,0]
a5 = [x1-la, x2-la, x3-la, 0,0,0,(xaI-la)*sin(theta),0,0,0,0,0]         #Marek Edit:
a6 = [0, 0, 0, 0,0,0,0, x1, 1,0,0,hh]
a7 = [0,0,0, 0, 0, 0, 0, 0, 0, x1, 1, 0]
a8 = [F1_Ry123(x2,x1), 0, 0, 0,0,0,F1_RaIaII(x2, xaI), x2, 1,0,0,hh]
a9 = [0,0,0, F2_Rz123(x2,x1), 0, 0, F2_RaIaII(x2, xaI), 0, 0, x2, 1, 0]
a10 = [F1_Ry123(x3,x1), F1_Ry123(x3,x2), 0, 0,0,0,F1_RaIaII(x3, xaI), x3, 1,0,0,hh]
a11 = [0,0,0, F2_Rz123(x3,x1), F2_Rz123(x3,x2), 0, F2_RaIaII(x3, xaI), 0, 0, x3, 1, 0]
a12 = [F3A_Ry123(xaI, x1)+F3B_Ry123(xaI, x1), 0, 0,
       F3B_Rz123(xaI, x1), 0, 0, 0, xaI*sin(theta), sin(theta), xaI*cos(theta),cos(theta),zh*sin(theta)-h/2*cos(theta)]

 
b1 = [-P*sin(theta) - Q1(la)]
b2 = [-P*cos(theta)]
b3 = [-P*(sin(theta)*zh-cos(theta)*h/2) - al.integration(1,la,al.qT_tilde,al.x_coor) ]      #Marek edit: replaced + by -
b4 = [P*cos(theta)*(la-xaII)]
b5 = [P*sin(theta)*(la-xaII) + al.integration(2,la,mult*al.q_tilde,al.x_coor)]
b6 = [-d1*cos(theta) - F1_const(x1) - 0]
b7 = [d1*sin(theta) - F2_const(x1)]
b8 = [- F1_const(x2) - 0]
b9 = [-0]
b10 = [-d3*cos(theta) - F1_const(x3)- F1_RaIaII(x3, xaII)*P]
b11 = [d3*sin(theta) - F2_const(x3)]
b12 = [-F3A_const(xaI) - F3B_const(xaI)-0]    #Marek Edit, added F3P

print("Hellooo", F3P_const(xaI,xaII))

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

### matrix equations verification
sumFy = Ry1 + Ry2 + Ry3 + RaI*sin(theta) + P*sin(theta) + Q1(la)
sumFz = Rz1 + Rz2 + Rz3 + RaI*cos(theta) + P*cos(theta)
sumMx = RaI*h*(sin(theta)-cos(theta))/2 - P*h*(cos(theta)-sin(theta))/2 - aero_Mx(la)
sumMy = (x1-x2)*Rz2 + (x1-x3)*Rz3 + (x1-xaI)*cos(theta)*RaI - P*cos(theta)*(xaII-x1)
sumMz = Ry1*(xaI-x1) + Ry2*(x2-xaI) + Ry3*(x3-xaI) + P*sin(theta)*(xaII-xaI) + aero_Mz(la)

BC1 = 1/(E*Izz)*Q4(x1) + hh/(G*J)*tao2(x1) + hh*C5 + C1*x1 + C2 + d1*cos(theta)
BC2 = C3*x1 + C4 -d1*sin(theta)
BC3 = F1_Ry123(x2,x1)*Ry1 + F1_RaIaII(x2, xaI)*RaI + 1/(E*Izz)*Q4(x2) + hh/(G*J)*tao2(x2) + C1*x2 + C2 + hh*C5
BC4 = 1/(6*E*Iyy)*( Rz1*Mac(x2-x1)**3 + RaI*cos(theta)*Mac(x2-xaI)**3 ) + C3*x2 + C4

BC5 = F1_Ry123(x3,x1)*Ry1 + F1_Ry123(x3,x2)*Ry2 + F1_RaIaII(x3,xaI)*RaI + F1_const(x3) + C1*x3 + C2 + hh*C5 + d3*cos(theta) + F1_RaIaII(x3, xaII)*P

BC6 = F2_Rz123(x3,x1)*Rz1 + F2_Rz123(x3,x2)*Rz2 + F2_RaIaII(x3,xaI)*RaI + F2_RaIaII(x3,xaII)*P  +C3*x3 + C4 - d3*sin(theta)

BC7 = F3A_Ry123(xaI, x1)*Ry1 + F3B_Ry123(xaI,x1)*Ry1 + F3B_Rz123(xaI,x1)*Rz1 + F3A_const(xaI) + F3B_const(xaI) + C1*xaI*sin(theta) + C2*sin(theta) + C3*xaI*cos(theta) + C4*cos(theta)+ (zh*sin(theta)-h/2*cos(theta))*C5
###


###
def v(x):
    return 1/(E*Izz)* ( 1/6 *Ry1*Mac(x-x1)**3 + 1/6 *Ry2*Mac(x-x2)**3 + 1/6 *Ry3*Mac(x-x3)**3 
              + 1/6 *RaI*sin(theta)*Mac(x-xaI)**3 + 1/6 *P*sin(theta)*Mac(x-xaII)**3  + Q4(x) ) + C1*x + C2
              
def w(x):
    return 1/(E*Iyy) * ( 1/6*Rz1*Mac(x-x1)**3 + 1/6*Rz2*Mac(x-x2)**3 + 1/6*Rz3*Mac(x-x3)**3 
              + 1/6*RaI*cos(theta)*Mac(x-xaI)**3 + 1/6*P*cos(theta)*Mac(x-xaII)**3 )  + C3*x + C4
    
def phi(x):
    return 1/(G*J) * ( Ry1*hh*Mac(x-x1) + Ry2*hh*Mac(x-x2) + Ry3*hh*Mac(x-x3) 
              + RaI*Mac(x-xaI)*(zh*sin(theta) - h/2*cos(theta)) 
              + P*Mac(x-xaII)*(zh*sin(theta) - h/2*cos(theta)) 
              + tao2(x) )  +  C5
 
import matplotlib.pyplot as plt

xx = np.linspace(0,la, 100)
plt.plot(xx,v(xx))
plt.show()

plt.plot(xx,w(xx))
plt.show()

plt.plot(xx,phi(xx))
plt.show()




dY = v(xx) + hh*phi(xx)

dYe = dY*cos(theta) - w(xx)*sin(theta)
dZe = w(xx)*cos(theta)+dY*sin(theta)



hinge = [x1,x2,x3]
aIaII = [xaI, xaII]
bcY = [d1,0,d3]
bcZ = [0,0,0]

#
def powzero(x):
    if x ==0:
        return 0
    else:
        return 1.
def Mz(x):
    M = np.zeros(np.shape(x))
    for i in range(len(x)):
        M[i] = -Ry1*Mac(x[i]-x1) - Ry2*Mac(x[i]-x2) - Ry3*Mac(x[i]-x3) - RaI*sin(theta)*Mac(x[i]-xaI) - P*sin(theta)*Mac(x[i]-xaII) - al.integration(2,x[i],mult*al.q_tilde,al.x_coor)
    return M

def Vy(x):
    V = np.zeros(np.shape(x))
    for i in range(len(x)):
        V[i] = -Ry1*powzero(Mac(x[i]-x1)) - Ry2*powzero(Mac(x[i]-x2)) - Ry3*powzero(Mac(x[i]-x3)) - RaI*sin(theta)*powzero(Mac(x[i]-xaI)) - P*sin(theta)*powzero(Mac(x[i]-xaII)) - al.integration(1,x[i],mult*al.q_tilde,al.x_coor)
    return V


def My(x):
    M = np.zeros(np.shape(x))
    for i in range(len(x)):
        M[i] = -Rz1*Mac(x[i]-x1) - Rz2*Mac(x[i]-x2) - Rz3*Mac(x[i]-x3) - RaI*cos(theta)*Mac(x[i]-xaI) - P*cos(theta)*Mac(x[i]-xaII)
    return M

def Vz(x):
    V = np.zeros(np.shape(x))
    for i in range(len(x)):
        V[i] = -Rz1*powzero(Mac(x[i]-x1)) - Rz2*powzero(Mac(x[i]-x2)) - Rz3*powzero(Mac(x[i]-x3)) - RaI*cos(theta)*powzero(Mac(x[i]-xaI)) - P*cos(theta)*powzero(Mac(x[i]-xaII))
    return V



plt.plot(xx,Vz(xx))
#
plt.plot(xx,Vy(xx))
plt.show()
#plt.plot(xx,phi(xx))


plt.plot(xx,dZe)
plt.scatter(hinge, bcZ)
plt.plot(xx,-dYe)
plt.scatter(hinge, bcY)
plt.scatter(aIaII,[0,0])
plt.show()


    
