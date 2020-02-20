import numpy as np
from variables import *
from booms_locations import *
from MoI import *

#new constants and knowns
ha = h/2                                    #Ha is the half of the aileron height!!
l_sk = np.sqrt((Ca-(ha))**2)+((ha)**2)  #Length of region 3 and 4
A_I = (1/2)*np.pi*(ha)**2                 #Area of cell I
A_II = (Ca-(ha))*(ha)                   #Area of cell II
qb01 = 0 #Due to cut
qb02 = 0 #Due to cut
qb05 = 0 #Due to cut


#Formulas base shear flows after analytical integration by hand
qb1 = -(1/Izz_Aileron)*(t_sk*ha**2+boom_area*(y_b1+y_b2))+qb01

qb2 = -(1/Izz_Aileron)*((1/2)*t_sp*ha**2)+qb02

qb03 = qb1+qb2

qb3 = -(1/Izz_Aileron)*(t_sk*ha*l_sk-(1/2)*(ha/l_sk)*l_sk**2+boom_area*(y_b3+y_b4+y_b5+y_b6))+qb03

qb04 = qb3
                
qb4 = -(1/Izz_Aileron)*(-(1/2)*t_sk*(ha/l_sk)*l_sk**2+boom_area*(y_b7+y_b8+y_b9+y_b10))+qb04

qb5 = -(1/Izz_Aileron)*((1/2)*t_sp*ha**2)+qb05

qb06 = qb4-qb5
                
qb6 = -(1/Izz_Aileron)*(-t_sk*ha**2+boom_area*y_b11)+qb06

print('qb1 =',qb1, 'qb2 =', qb2, 'qb3 =', qb3,'qb4 =', qb4,'qb5 =', qb5,'qb6 =', qb6)

# Shear flows from the two cells
# Solve Matrix for qsoI and qsoII

A = np.array([[(np.pi*h/t_sk)+(2*h/t_sp), -(2*h/t_sp)],
             [-2*h/t_sp, (2*h/t_sp)+(2*l_sk/t_sk)]])
B = np.array([(-qb1*h*np.pi)/(2*t_sk)-(qb6*h*np.pi)/(2*t_sk)-(qb2*h)/(t_sp)+(qb5*h)/(t_sp),
             (-qb2*h/t_sp)-(qb5*h/t_sp)-(qb3*l_sk/t_sk)+(qb4*l_sk/t_sk)])
M = np.linalg.solve(A, B)

print('qbsoI =', M[0], 'qsoII =', M[1])


