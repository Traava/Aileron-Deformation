import numpy as np
from variables import * 
from math import *
import matplotlib.pyplot as plt



l_skin           = sqrt((0.5*h)**2 + (Ca-0.5*h)**2)          #length of the inclined side of one side aileron, due to symmetry the length of the top and bottom part are the same
theta_skin       = atan((0.5*h)/(Ca-0.5*h))                  #angle of the inclined sides with the horizontal
boom_area        = (h_st+w_st)*t_st
booms_separation = (2*l_skin + pi*(0.5*h))/11                #circumfrance of the entire structure summed and divided by number of stiffeners

#y-coordinates of all the booms, with the z-axis through the symmetry plan of the corss-section
#booms are numbered from the one in the leading edge and clockwise
theta_separation = booms_separation/(h/2)                    #angle of the second boom  with the z-axis, measured from the center of the semi-circular part
y_b1  = 0
y_b2  = (0.5*h)*sin(theta_separation)
y_b3  = 0.5*h - (2*booms_separation - 0.5*pi*(0.5*h))*sin(theta_skin)
y_b4  = y_b3 - booms_separation*sin(theta_skin)
y_b5  = y_b4 - booms_separation*sin(theta_skin)
y_b6  = y_b5 - booms_separation*sin(theta_skin)
y_b7  = -y_b6  # Because of symmetry
y_b8  = -y_b5  # Because of symmetry
y_b9  = -y_b4  # Because of symmetry
y_b10 = -y_b3  # Because of symmetry
y_b11 = -y_b2  # Because of symmetry

y_coordinates = np.array([y_b1,y_b2,y_b3,y_b4,y_b5,y_b6,y_b7,y_b8,y_b9,y_b10,y_b11])    #array containing all the y-coordinates


#z-coordinates of the booms measured from the leading edge of the aileron
#this values are used to plot the aileron with the booms, but they are adjusted later to the values with respect to the centroid
z_b1  = 0
z_b2  = (h/2)*(1-cos(theta_separation))
z_b3  = (h/2) + (2*booms_separation - 0.5*pi*(0.5*h))*cos(theta_skin)
z_b4  = z_b3 + booms_separation*cos(theta_skin)
z_b5  = z_b4 + booms_separation*cos(theta_skin)
z_b6  = z_b5 + booms_separation*cos(theta_skin)
z_b7  = z_b6
z_b8  = z_b5
z_b9  = z_b4
z_b10 = z_b3
z_b11 = z_b2

z_coordinates = np.array([z_b1,z_b2,z_b3,z_b4,z_b5,z_b6,z_b7,z_b8,z_b9,z_b10,z_b11])    #array containing all the z-coordinates

    

#circle equation
z1 = np.linspace(0, 0.5*h, 50)
y1 = np.sqrt((0.5*h)**2 - (z1-0.5*h)**2)      #top quarter
y11 = -y1                                     #bottom quarter


#lines equation
slope     = (-0.5*h)/(Ca-(0.5*h))
intercept = -slope*Ca
z2 = np.linspace(0.5*h, Ca, 150)
y2 = slope*z2 + intercept                   #top skin
y22 = -y2                                   #lower skin


if __name__ == '__main__': #only plot if running this file directly
    plt.plot(z1,y1, 'k')                        #plotting top quarter
    plt.plot(z1,y11, 'k')                       #plotting bottom quarter
    plt.plot(z2,y2, 'k')                        #plotting top skin
    plt.plot(z2,y22, 'k')                       #plotting bottom skin
    plt.vlines(x=0.5*h, ymin=-0.5*h, ymax=0.5*h, linewidth = 3)    #plotting spar 
    plt.plot(z_coordinates,y_coordinates,'ro')   #plotting Booms



#Centroid Calculations
sum = 0 

for i in range(1,12):             #multiplying each boom by it's z-coordinate
    j = eval('z_b'+ str(i))
    sum = sum + j*boom_area 

centroid_semicircle = ((h/2) - (h/pi))*(pi*0.5*h*t_sk)       #centroid of semi-circular part multiplied by its area
centroid_skin       = (0.5*Ca + 0.25*h)*(t_sk*l_skin)        #centroid of skin multiplied by its area
centroid_spar       = (0.5*h)*(h*t_sp)                       #centroid of the spar part multiplied by its area

Area = boom_area*11 + (h*t_sp) + (t_sk*l_skin*2) + (pi*0.5*h*t_sk)      #total area of the cross-section

centroid = (sum + 2*centroid_skin + centroid_semicircle + centroid_spar)/Area   #Final z-coordinate of the centroid with respect to the leading edge


if __name__ == '__main__':  
    plt.scatter(x = centroid, y=0)                     #plotting the centroid location
    plt.show()



#Adjusted z-coordinates of the cross-section with respect to the centroid not the leading edge
#these coordinated will be used in calculating the contribution of booms in the shear flow
z_coordinates_adjusted = np.array([])

for i in range(1,12):
    j = eval('z_b'+str(i))
    z_coor = j - centroid
    z_coordinates_adjusted = np.append(z_coordinates_adjusted, z_coor)
    
