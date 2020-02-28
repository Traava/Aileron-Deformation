import numpy as np
from variables import * 
from math import *
import matplotlib.pyplot as plt
from MatrixSystem import *
from MoI import *

    
n    = 300   #half the Number of panels the circumference is divided into, needs to be an even number
n_sp = 100
m    = 80    #m is the number of sections the span is divided into

l_skin = sqrt((0.5*h)**2 + (Ca-0.5*h)**2)   #length of the inclined side of one side aileron, due to symmetry the length of the top and bottom part are the same
panel_sep = (2*l_skin + pi*(0.5*h))/(2*n)   #Distance between the panels

y_coord = np.array([])
z_coord = np.array([])

## still add panel coordinates opposite side
#coordinates with respect to the leading edge for half of the circumference
for i in range(n):
    if i*panel_sep <= h*pi/4:                       #panel is on the arc part of the circumference
        theta_sep = i*panel_sep/(h/2)               #angle of the second boom  with the z-axis, measured from the center of the semi-circular part
        zp_i = (h/2)*(1-cos(theta_sep)) - centroid
        yp_i = -(0.5*h)*sin(theta_sep)
        y_coord = np.append(y_coord,yp_i)
        y_coord = np.append(y_coord,-yp_i)
        z_coord = np.append(z_coord,zp_i)
        z_coord = np.append(z_coord,zp_i)
        
    elif i*panel_sep > h*pi/4:                      #panel is on the triangular part of the cross section
        theta_skin = atan((0.5*h)/(Ca-0.5*h))       #angle of the inclined sides with the horizontal
        zp_i = ((h/2) + (i*panel_sep - 0.5*pi*(0.5*h))*cos(theta_skin)) - centroid
        yp_i = - (0.5*h - (i*panel_sep - 0.5*pi*(0.5*h))*sin(theta_skin))
        y_coord = np.append(y_coord,yp_i)
        y_coord = np.append(y_coord,-yp_i)
        z_coord = np.append(z_coord,zp_i)
        z_coord = np.append(z_coord,zp_i)


########SPAR########
z_coord_sp = (h/2 - centroid)*np.ones(n_sp)
y_coord_sp = np.zeros(n_sp)
sep = h/n_sp
for i in range(n_sp):
        s = -h/2 + sep*i
        y_coord_sp[i] = s



#Bending stresses
#from MatrixSystem import Mz, My
#Dividing the span into sections and calculating the maximum moment at each section.

x_sectionstart = np.array([0])
for j in range(1,m,1): #defining start of each section
    x_sectionstart = np.append(x_sectionstart,j*(la/m))  #la is span of aileron, defined in variables
    #j = j+1
x_sectionend =  np.array([la/m])
for k in range(2,m+1,1):#defining end of each section
    x_sectionend = np.append(x_sectionend,(k+1)*(la/m)) 
    #k=k+1
#Calculate the maximum moment at each section as this will give maximum/critical stress
#Mz_section = np.array([]) #making an empty array to which the maximum moment at each section will be appended
#My_section = np.array([]) #making an empty array to which the maximum moment at each section will be appended

#for l in range(0,m):
    
steps = (x_sectionstart + x_sectionend)/2

M_z = Mz(steps)
M_y = My(steps)

Stress = np.zeros((m,2*n))
for i in range(m):
    for j in range(2*n):
        sigma = M_z[i]*y_coord[j]/Izz_Aileron + M_y[i]*z_coord[j]/Iyy_Aileron
        Stress[i][j] = sigma


Stress_sp = np.zeros((m,n_sp))
for i in range(m):
    for j in range(n_sp):
        sigma_sp = M_z[i]*y_coord_sp[j]/Izz_Aileron + M_y[i]*z_coord_sp[j]/Iyy_Aileron
        Stress_sp[i][j] = sigma_sp


print('x =', (la/m)*(20-1)+0.5*sep)
print(np.max(Stress))
print(np.max(Stress_sp))

Max_location = np.where(Stress == np.max(np.amax(Stress, axis = 0)))
Max_loc = np.asscalar(Max_location[0])
#loc = (Max_loc-1)*la/n

print(Max_location)
print(Max_loc)



y_coord = np.hstack((y_coord,y_coord_sp))
z_coord = np.hstack((z_coord,z_coord_sp))
Stress = np.hstack((Stress[Max_loc],Stress_sp[Max_loc]))


fig1 = plt.figure()
ax1 = plt.axes()
ax1.invert_yaxis()
img1 = plt.scatter(z_coord, y_coord, c = Stress, cmap ='RdYlGn_r')
#img2 = plt.scatter(z_coord_sp, y_coord_sp, c = Stress_sp[20], cmap ='RdYlGn_r') 


cbar1 = fig1.colorbar(img1)
cbar1.set_label(r'[$\frac{N}{m^2}$]')
ax1.set_title("Direct Stress Distribution")
ax1.set_xlabel('z [m]')
ax1.set_ylabel('y [m]')

