import numpy as np
from variables import * 
from math import *
import matplotlib.pyplot as plt

    
n = 20                                      #half the Number of panels the circumference is divided into, needs to be an even number
l_skin = sqrt((0.5*h)**2 + (Ca-0.5*h)**2)   #length of the inclined side of one side aileron, due to symmetry the length of the top and bottom part are the same
panel_sep = (2*l_skin + pi*(0.5*h))/(2*n)   #Distance between the panels

y_coord = np.array(0)
z_coord = np.array(0)

## still add panel coordinates opposite side
#coordinates with respect to the leading edge for half of the circumference
for i in range(1,n+1):
    if i*panel_sep <= h*pi/4:                       #panel is on the arc part of the circumference
        theta_sep = i*panel_sep/(h/2)               #angle of the second boom  with the z-axis, measured from the center of the semi-circular part
        zp_i = (h/2)*(1-cos(theta_sep))
        yp_i = (0.5*h)*sin(theta_sep)
        np.append(y_coord,yp_i)
        np.append(z_coord,zp_i)
        i = i+1
    elif i*panel_sep > h*pi/4:                      #panel is on the triangular part of the cross section
        theta_skin = atan((0.5*h)/(Ca-0.5*h))       #angle of the inclined sides with the horizontal
        zp_i = (h/2) + (i*panel_sep - 0.5*pi*(0.5*h))*cos(theta_skin)
        yp_i = 0.5*h - (i*panel_sep - 0.5*pi*(0.5*h))*sin(theta_skin)
        np.append(y_coord,yp_i)
        np.append(z_coord,zp_i)
        i = i+1

#circle equation
z1 = np.linspace(0, 0.5*h, 50)
y1 = np.sqrt((0.5*h)**2 - (z1-0.5*h)**2)      #top quarter
y11 = -y1                                     #bottom quarter


#lines equation
slope     = (-0.5*h)/(Ca-(0.5*h))
intercept = -slope*Ca
z2 = np.linspace(0.5*h, Ca, 150)
y2 = slope*z2 + intercept                                       #top skin
y22 = -y2                                                       #lower skin

plt.plot(z1,y1, 'k')                                            #plotting top quarter
plt.plot(z1,y11, 'k')                                           #plotting bottom quarter
plt.plot(z2,y2, 'k')                                            #plotting top skin
plt.plot(z2,y22, 'k')                                           #plotting bottom skin
plt.vlines(x=0.5*h, ymin=-0.5*h, ymax=0.5*h, linewidth = 3)     #plotting spar 
plt.plot(z_coord,y_coord,'ro')                                  #plotting panels
plt.show()


#Bending stresses
#from MatrixSystem import Mz(x), My(x)
#la is span of aileron, defined in variables
#k= np.arange(0,m,1) # where k gives the section number and m the number of sections
#xsectionlst = np.array()
#for i in k: # creating an array with coordinates of each section
#    xsection = (la/m)*(i+1/2) # xcoordinates at which the moments of each cross section are evaluated are in the middle of each section
#    np.append(xsectionlst, xsection)
#from MoI import Izz_Aileron, Iyy_Aileron
#sigmax = np.array(0)
#bending_stresses = np.array(0)
#for x in xsectionlst:          #k is number of sections along the span
#    for i in range(1,n+1):      #i is each panel along the cross section
#        sigmax_i = Mz[x]*y_coord[i]/Izz_Aileron + My[x] z_coord[i]/Iyy_Aileron
#        np.append(sigmax,sigmax_i)
#        i= i+1
#    np.append(sigmax,bending_stresses)
#    k=k+1








