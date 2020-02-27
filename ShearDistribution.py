from math import *
import numpy as np
import matplotlib.pyplot as plt
from MoI import Izz_Aileron, Iyy_Aileron
from variables import *
from booms_locations import *

#Numerical integrator that gives an array that calculates the values of the integral
#from the initial point till the different point on the contour 

def integrate(fn, x):
    dx = (x[-1]-x[0])/len(x)
    
    integral = np.array([])
    f0 = fn[0]

    for i in range(len(x)):
        
        integrate_temp = 0
        integrate = 0
        f_n = fn[i]
        
        for j in range(1, i):
            integrate_temp += fn[j]
            
        integrate += (dx/2)*(f0 + f_n + 2*integrate_temp)
        integral = np.append(integral, integrate )
        
    return(integral)
    
    
r = h/2
l_sk = sqrt((0.5*h)**2 + (Ca-0.5*h)**2)
A_I  = (pi*r**2)/2
A_II = r*(Ca-r)
     



##---------------------Vertical Shear flow calculations----------------------#
def VerticalShear(Sy):
    qb01 = 0
    qb02 = 0
    qb05 = 0    

##_--------------------------region 1----------------------------------#
    def function(x):
        return np.sin(x)
    
    x1 = np.linspace(0,theta_separation,50)
    qb11 = -(1/Izz_Aileron)*((-(integrate(function(x1),x1))*r**2*t_sk)) +qb01
    x2 = np.linspace(theta_separation,pi/2,50)
    qb12 = -(1/Izz_Aileron)*((-(integrate(function(x2),x2))*r**2*t_sk) + boom_area*(y_b1+y_b2)) + qb11[-1]
    
    qb1 = np.hstack((qb11, qb12))
    
    
    ##_--------------------------region 2----------------------------------#
    def function(x):
        return x
    x = np.linspace(0,-r,50)
    qb2 = -(1/Izz_Aileron)*((integrate(function(x),x))*t_sp) + qb02
    
    
    ##_--------------------------region 3----------------------------------#
    def function(x):
        return -1 + x/l_sk
    
    separation_initial = 2*booms_separation-pi*r/2
    x = np.linspace(0,separation_initial,50)
    qb03 = qb1[-1] + qb2[-1]
    qb31 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk)) + qb03
    
    x = np.linspace(separation_initial,booms_separation+separation_initial,50)
    qb32 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(y_b3)) + qb31[-1]
    
    x = np.linspace(booms_separation+separation_initial,2*booms_separation+separation_initial,50)
    qb33 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(y_b4)) + qb32[-1]
    
    x = np.linspace(2*booms_separation+separation_initial,3*booms_separation+separation_initial,50)
    qb34 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(y_b5)) + qb33[-1]
    
    x = np.linspace(3*booms_separation+separation_initial,l_sk,50)
    qb35 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(y_b6)) + qb34[-1]
    
    qb3 = np.hstack((qb31,qb32,qb33,qb34,qb35))
    
    
    
    ##_--------------------------region 4----------------------------------#
    def function(x):
        return x
    #x = np.linspace(0,l_sk,50)
    qb04 = qb3[-1]
    #qb4 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(y_b7+y_b8+y_b9+y_b10)) + qb04
    
    
    separation_initial = booms_separation/2
    
    x = np.linspace(0,separation_initial,50)
    qb41 = -(1/Izz_Aileron)*(integrate(function(x),x)*r*t_sk/l_sk) + qb04
    
    x = np.linspace(separation_initial,booms_separation+separation_initial,50)
    qb42 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(y_b7)) + qb41[-1]
    
    x = np.linspace(booms_separation+separation_initial,2*booms_separation+separation_initial,50)
    qb43 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(y_b8)) + qb42[-1]
    
    x = np.linspace(2*booms_separation+separation_initial,3*booms_separation+separation_initial,50)
    qb44 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(y_b9)) + qb43[-1]
    
    x = np.linspace(3*booms_separation+separation_initial,l_sk,50)
    qb45 = -(1/Izz_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(y_b10)) + qb44[-1]
    
    qb4 = np.hstack((qb41,qb42,qb43,qb44,qb45))
    
    
    ##_--------------------------region 5----------------------------------#
    def function(x):
        return x
    x = np.linspace(0,r,50)
    qb5 = -(1/Izz_Aileron)*((integrate(function(x),x))*t_sp) + qb05
    
    
    ##_--------------------------region 6----------------------------------#
    def function(x):
        return np.sin(x)
    qb06 = qb4[-1] + qb5[-1]
    
    x = np.linspace(-pi/2,-theta_separation,50)
    qb61 = -(1/Izz_Aileron)*((integrate(function(x),x))*r**2*t_sk)  + qb06
    x = np.linspace(-theta_separation,0,50)
    qb62 = -(1/Izz_Aileron)*((integrate(function(x),x))*r**2*t_sk + boom_area*(y_b11))  + qb61[-1]
    
    qb6 = np.hstack((qb61,qb62))
    
    
    a11 = (np.pi*r/t_sk)+(2*r/t_sp)
    a12 = -(2*r/t_sp)
    a21 = -2*r/t_sp
    a22 = (2*r/t_sp)+(2*l_sk/t_sk)
    
    A = np.array([[a11, a12],
                 [a21, a22]])
        
    b11 = ((-(qb1[-1]+qb6[-1])*np.pi*r)/(2*t_sk))- (((qb5[-1]-qb2[-1])*r)/t_sp)
    b21 = ((-(qb3[-1]+qb4[-1])*l_sk)/t_sk) - (((qb2[-1]-qb5[-1])*r)/t_sp)
    B = np.array([[b11],
                   [b21]])
    M = np.linalg.solve(A, B)
    
    q1 = Sy*(qb1+M[0])
    q2 = Sy*(qb2-M[0]+M[1])
    q3 = Sy*(qb3+M[1])
    q4 = Sy*(qb4+M[1])
    q5 = Sy*(qb5+M[0]-M[1])
    q6 = Sy*(qb6+M[0])
    
    #shear_center = (q1[-1]*np.pi*ha**2/2)+(q6[-1]*np.pi*ha**2/2)+(q3[-1]*l_sk*M_arm34)+(q4[-1]*l_sk*M_arm34)
    #zsc = np.asscalar(shear_center) + ha
    
    return q1, q2, q3, q4, q5, q6


##---------------------------------------------------------------------------------------------#
##---------------------------------------------------------------------------------------------#
##---------------------------------------------------------------------------------------------#
##-----------------------Horizontal Shear flow calculations------------------------------------#
def HorizontalShear(Sz):
    
    qbz01 = 0                                   #Due to cut
    qbz02 = 0                                   #Due to cut
    qbz05 = 0                                   #Due to cut
    
    
    ##_--------------------------region 1----------------------------------#
    def function(x):
        return (r*(1-np.cos(x))-centroid)
    
    x1 = np.linspace(0,theta_separation,50)
    qb11 = -(1/Iyy_Aileron)*((integrate(function(x1),x1))*r*t_sk) + qbz01
    x2 = np.linspace(theta_separation,pi/2,50)
    qb12 = -(1/Iyy_Aileron)*(((integrate(function(x2),x2))*r*t_sk) + boom_area*(z_b1+z_b2)) + qb11[-1]
    
    qbz1 = np.hstack((qb11, qb12))
    
    ##_--------------------------region 2----------------------------------#
    x = np.linspace(0,-r,50)
    y = np.ones(len(x))*(r-centroid)
    qbz2 = -(1/Iyy_Aileron)*((integrate(y,x))*t_sp) + qbz02
    
    ##_--------------------------region 3----------------------------------#
    def function(x):
        return (-(centroid-r) + (Ca-r)*x/l_sk)
    
    separation_initial = 2*booms_separation-pi*r/2
    x = np.linspace(0,separation_initial,50)
    qbz03 = qbz1[-1] + qbz2[-1]
    qb31 = -(1/Iyy_Aileron)*((integrate(function(x),x))*r*t_sk) + qbz03
    
    x = np.linspace(separation_initial,booms_separation+separation_initial,50)
    qb32 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(z_b3)) + qb31[-1]
    
    x = np.linspace(booms_separation+separation_initial,2*booms_separation+separation_initial,50)
    qb33 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(z_b4)) + qb32[-1]
    
    x = np.linspace(2*booms_separation+separation_initial,3*booms_separation+separation_initial,50)
    qb34 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(z_b5)) + qb33[-1]
    
    x = np.linspace(3*booms_separation+separation_initial,l_sk,50)
    qb35 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk) + boom_area*(z_b6)) + qb34[-1]
    
    qbz3 = np.hstack((qb31,qb32,qb33,qb34,qb35))
    
    
    ##_--------------------------region 4----------------------------------#
    def function(x):
        return ((Ca-centroid) - (Ca-r)*x/l_sk)

    qb04 = qbz3[-1]  
  
    separation_initial = booms_separation/2
    
    x = np.linspace(0,separation_initial,50)
    qb41 = -(1/Iyy_Aileron)*(integrate(function(x),x)*r*t_sk/l_sk) + qb04
    
    x = np.linspace(separation_initial,booms_separation+separation_initial,50)
    qb42 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(z_b7)) + qb41[-1]
    
    x = np.linspace(booms_separation+separation_initial,2*booms_separation+separation_initial,50)
    qb43 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(z_b8)) + qb42[-1]
    
    x = np.linspace(2*booms_separation+separation_initial,3*booms_separation+separation_initial,50)
    qb44 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(z_b9)) + qb43[-1]
    
    x = np.linspace(3*booms_separation+separation_initial,l_sk,50)
    qb45 = -(1/Iyy_Aileron)*(((integrate(function(x),x))*r*t_sk/l_sk) + boom_area*(z_b10)) + qb44[-1]
    
    qbz4 = np.hstack((qb41,qb42,qb43,qb44,qb45))

    ##_--------------------------region 5----------------------------------#
    x = np.linspace(0,r,50)
    y = np.ones(len(x))*(r-centroid)
    qbz5 = -(1/Iyy_Aileron)*((integrate(function(x),x))*t_sp) + qbz05
    
    ##_--------------------------region 6----------------------------------#
    def function(x):
        return (r*(1-np.cos(x)) - centroid)
    qb06 = qbz4[-1] + qbz5[-1]
    
    x = np.linspace(-pi/2,-theta_separation,50)
    qb61 = -(1/Iyy_Aileron)*((integrate(function(x),x))*r*t_sk)  + qb06
    x = np.linspace(-theta_separation,0,50)
    qb62 = -(1/Iyy_Aileron)*((integrate(function(x),x))*r*t_sk + boom_area*(z_b11))  + qb61[-1]
    
    qbz6 = np.hstack((qb61,qb62))
    
    
    
    qz1 = Sz*qbz1
    qz2 = Sz*qbz2
    qz3 = Sz*qbz3
    qz4 = Sz*qbz4
    qz5 = Sz*qbz5
    qz6 = Sz*qbz6
    
    return qz1, qz2, qz3, qz4, qz5, qz6


def Torsion(T):
    
    b= np.array([[0],[0],[T]])


    a11 = (h/(2*A_I))*((pi/(2*t_sk))+(1/t_sp))
    a12 = -h/(2*A_I*t_sp)
    a13 = -1
    a21 = -h/(2*A_II*t_sp)
    a22 = (1/(2*A_II))*((2*l_skin/t_sk)+(h/t_sp))
    a23 = -1
    a31 = 2*A_I
    a32 = 2*A_II
    a33 = 0
    
    A = np.array([[a11, a12, a13],
                  [a21, a22, a23],
                  [a31, a32, a33]])
  
    x = np.linalg.solve(A,b)
    
    
    #J = np.asscalar(T/x[2]) #Torsional stiffness calculation, used only for verification purposes
    
    qb1 = np.ones(100)*x[0]
    qb2 = np.ones(50)*(x[1]-x[0])
    qb3 = np.ones(250)*x[1]
    qb4 = np.ones(250)*x[1]
    qb5 = np.ones(50)*(x[0]-x[1])
    qb6 = np.ones(100)*x[0]
    
    
    return qb1,qb2,qb3,qb4,qb5,qb6



def ShearDistribution(Torque, Vertical_shear, Horizontal_shear):
    qb1 = VerticalShear(Vertical_shear)[0] + HorizontalShear(Horizontal_shear)[0] + Torsion(Torque)[0]
    qb2 = VerticalShear(Vertical_shear)[1] + HorizontalShear(Horizontal_shear)[1] + Torsion(Torque)[1]
    qb3 = VerticalShear(Vertical_shear)[2] + HorizontalShear(Horizontal_shear)[2] + Torsion(Torque)[2]
    qb4 = VerticalShear(Vertical_shear)[3] + HorizontalShear(Horizontal_shear)[3] + Torsion(Torque)[3]
    qb5 = VerticalShear(Vertical_shear)[4] + HorizontalShear(Horizontal_shear)[4] + Torsion(Torque)[4]
    qb6 = VerticalShear(Vertical_shear)[5] + HorizontalShear(Horizontal_shear)[5] + Torsion(Torque)[5]
    
    return qb1,qb2,qb3,qb4,qb5,qb6




























