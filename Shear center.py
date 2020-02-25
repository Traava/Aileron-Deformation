import numpy as np
from variables import *
from booms_locations import *
from MoI import *

###new constants and knowns
ha = h/2                                    #Ha is the half of the aileron height
#l_sk = np.sqrt((Ca-(ha))**2)+((ha)**2)      #Length of region 3 and 4
l_sk = l_skin
A_I = (1/2)*np.pi*(ha)**2                   #Area of cell I
A_II = (Ca-(ha))*(h)                       #Area of cell II
qb01 = 0                                    #Due to cut
qb02 = 0                                    #Due to cut
qb05 = 0                                    #Due to cut


###Formulas base shear flows after analytical integration by hand
qb1 = -(1/Izz_Aileron)*(-t_sk*ha**2+boom_area*(y_b1+y_b2))+qb01
qb2 = -(1/Izz_Aileron)*((1/2)*t_sp*ha**2)+qb02
qb03 = qb1+qb2
qb3 = -(1/Izz_Aileron)*(-0.5*t_sk*ha*l_sk+boom_area*(y_b3+y_b4+y_b5+y_b6))+qb03
qb04 = qb3               
qb4 = -(1/Izz_Aileron)*(0.5*t_sk*l_sk*ha+boom_area*(y_b7+y_b8+y_b9+y_b10))+qb04
qb5 = -(1/Izz_Aileron)*((1/2)*t_sp*ha**2)+qb05
qb06 = qb4+qb5
qb6 = -(1/Izz_Aileron)*(-t_sk*ha**2+boom_area*y_b11)+qb06

#print('qb1 =',qb1, 'qb2 =', qb2, 'qb3 =', qb3,'qb4 =', qb4,'qb5 =', qb5,'qb6 =', qb6)


###Shear flows from the two cells: solve Matrix for qsoI and qsoII
a11 = (np.pi*ha/t_sk)+(2*ha/t_sp)
a12 = -(2*ha/t_sp)
a21 = -2*ha/t_sp
a22 = (2*ha/t_sp)+(2*l_sk/t_sk)

A = np.array([[a11, a12],
             [a21, a22]])
#B = np.array([(qb1*ha*np.pi)/(2*t_sk)+(qb6*ha*np.pi)/(2*t_sk)+(qb2*ha)/(t_sp)+(qb5*ha)/(t_sp),
#             (qb2*ha/t_sp)+(qb5*ha/t_sp)-(qb3*l_sk/t_sk)-(qb4*l_sk/t_sk)])

b11 = ((-(qb1+qb6)*np.pi*ha)/(2*t_sk))- (((qb5-qb2)*ha)/t_sp)
b21 = ((-(qb3+qb4)*l_sk)/t_sk) - (((qb2-qb5)*ha)/t_sp)
B = np.array([[b11],
               [b21]])
M = np.linalg.solve(A, B)

#print('qbsoI =', M[0], 'qsoII =', M[1])


###Final shear flows in each region due to unit shear force Sy
q1 = qb1+M[0]
q2 = qb2-M[0]+M[1]
q3 = qb3+M[1]
q4 = qb4+M[1]
q5 = qb5+M[0]-M[1]
q6 = qb6+M[0]

#print('q1 =',q1, 'q2 =', q2, 'q3 =', q3,'q4 =', q4,'q5 =', q5,'q6 =', q6)


###Calculate the shear center by setting Mi = Me = -shearcenter
###moment arm of regio 3 and 4
M_arm34 = ha*np.cos(np.arctan(ha/(Ca-ha))) 

###Shearcenter calculation
shear_center = (q1*np.pi*ha**2/2)+(q6*np.pi*ha**2/2)+(q3*l_sk*M_arm34)+(q4*l_sk*M_arm34)

zsc = np.asscalar(shear_center) + ha

###Since this shear center is measured from halfway the spar, and our shear center is
### measured from the leading edge, we have to add the radius of the semicircle!

print('Shear center distance from the leading edge =', zsc)

###Goede antwoord is 0.08553893540215983 van verification model
