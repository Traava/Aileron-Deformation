"This program will calculate the tosional stiffness J of the cross section."
"The cross section is divided into a semicircular left cell and a triangular right cell"

import numpy as np
from variables import * 
from math import *
from scipy import linalg, mat
from booms_locations import l_skin
#Step 1: Calculate the enclosed area's
A_1 =  (pi *(h/2)**2)/2 #The left enclosed area is a semicircle
A_2 = 0.5 * h * (Ca-(h/2))

#Step 2: Setting up equilibrium equations in a matrix and solve for q,01, q02 and Gdtheta/dz

b= np.array([[0],[0],[1]])

A = np.array([ [(1/(2*A_1))*(pi*h/(t_sk * 2)+h/t_sp),-h/(2*A_1*t_sp) , -1], [-h/(2*A_2*t_sp), (1/2*A_2)*((2*l_skin/t_sk )+h/t_sp),-1], [2*A_1, 2*A_2, 0] ]) 

x = np.linalg.solve(A,b)

##print("q_01 =",x[0]) 
##print("q_02 =",x[1]) 
##print("Gdtdz  =",x[2]) 


#Step 3: Calculate J with J= 1/Gdtheta/dz

print ("J =", 1/x[2])
