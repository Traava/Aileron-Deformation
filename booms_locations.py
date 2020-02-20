import numpy as np
from variables import * 
from math import *
import matplotlib.pyplot as plt

l_skin           = sqrt((0.5*h)**2 + (Ca-0.5*h)**2)          #length of the inclined side of the skin, top and bottom
theta_skin       = atan((0.5*h)/(Ca-0.5*h))                  #angle of the inclined sides with the horizontal
boom_area        = (h_st+w_st)*t_st
booms_separation = (2*l_skin + pi*(0.5*h))/11                #circumfrance of the entire structure summed and divided by number of stiffeners

#y-coordinates of all the booms, with the x-axis through the symmetry plan of the corss-section
#booms are numbered from the one in the leading edge and clockwise
theta_separation = booms_separation/(h/2)                    #angle of the second boom  with the x-axis, measured from the center of the semi-circular part
y_b1  = 0
y_b2  = (0.5*h)*sin(theta_separation)
y_b3  = 0.5*h - (2*booms_separation - 0.5*pi*(0.5*h))*sin(theta_skin)
y_b4  = y_b3 - booms_separation*sin(theta_skin)
y_b5  = y_b4 - booms_separation*sin(theta_skin)
y_b6  = y_b5 - booms_separation*sin(theta_skin)
y_b7  = -y_b6
y_b8  = -y_b5
y_b9  = -y_b4
y_b10 = -y_b3
y_b11 = -y_b2

y_coordinates = np.array([y_b1,y_b2,y_b3,y_b4,y_b5,y_b6,y_b7,y_b8,y_b9,y_b10,y_b11])    #array containing all the y-coordinates


#x-coordinates of the booms measured from the leading edge of the aileron
#this values are used to plot the aileron with the booms, but they are adjusted later to the values with respect to the centroid
x_b1  = 0
x_b2  = (h/2)*(1-cos(theta_separation))
x_b3  = (h/2) + (2*booms_separation - 0.5*pi*(0.5*h))*cos(theta_skin)
x_b4  = x_b3 + booms_separation*cos(theta_skin)
x_b5  = x_b4 + booms_separation*cos(theta_skin)
x_b6  = x_b5 + booms_separation*cos(theta_skin)
x_b7  = x_b6
x_b8  = x_b5
x_b9  = x_b4
x_b10 = x_b3
x_b11 = x_b2

x_coordinates = np.array([x_b1,x_b2,x_b3,x_b4,x_b5,x_b6,x_b7,x_b8,x_b9,x_b10,x_b11])    #array containing all the x-coordinates

#circle equation
x1 = np.linspace(0, 0.5*h, 50)
y1 = np.sqrt((0.5*h)**2 - (x1-0.5*h)**2)      #top quarter
y11 = -y1                                     #bottom quarter


#lines equation
slope     = (-0.5*h)/(Ca-(0.5*h))
intercept = -slope*Ca
x2 = np.linspace(0.5*h, Ca, 150)
y2 = slope*x2 + intercept                   #top skin
y22 = -y2                                   #lower skin

plt.plot(x1,y1, 'k')                        #plotting top quarter
plt.plot(x1,y11, 'k')                       #plotting bottom quarter
plt.plot(x2,y2, 'k')                        #plotting top skin
plt.plot(x2,y22, 'k')                       #plotting bottom skin
plt.vlines(x=0.5*h, ymin=-0.5*h, ymax=0.5*h, linewidth = 3)    #plotting spar 
plt.plot(x_coordinates,y_coordinates,'ro')   #plotting Booms



#Centroid Calculations
sum = 0 

for i in range(1,12):             #multiplying each boom by it's x-coordinate
    j = eval('x_b'+ str(i))
    sum = sum + j*boom_area 

centroid_semicircle = ((h/2) - (h/pi))*(pi*0.5*h*t_sk)       #centroid of semi-circular part multiplied by its area
centroid_skin       = (0.5*Ca + 0.25*h)*(t_sk*l_skin)        #centroid of skin multiplied by its area
centroid_spar       = (0.5*h)*(h*t_sp)                       #centroid of the spar part multiplied by its area

Area = boom_area*11 + (h*t_sp) + (t_sk*l_skin*2) + (pi*0.5*h*t_sk)      #total area of the cross-section

centroid = (sum + 2*centroid_skin + centroid_semicircle + centroid_spar)/Area   #Final x-coordinate of the centroid


plt.scatter(x = centroid, y=0)                     #plotting the centroid location

plt.show()



#Adjusted x-coordinates of the cross-section with respect to the centroid not the leading edge
#these coordinated will be used in calculating the contribution of booms in the shear flow
x_coordinates_adjusted = np.array([])

for i in range(1,12):
    j = eval('x_b'+str(i))
    x_coor = j - centroid
    x_coordinates_adjusted = np.append(x_coordinates_adjusted, x_coor)
    



