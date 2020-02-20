" Code calculating the moment of inertia of the cross section of the aileron"
" General approach: Dividing up the cross section into 4 parts; "
" the semicircular section, the spar, the inclined panel and the stiffeners idealized as booms"
" Then calculating the moi of each section around its centroid, then adding steiner terms wrt to the centroid of the whole cross section."

import numpy as np
from variables import * 
from math import *



#Step 1: Calculating the Z- and Y- Coordinates of the booms and centroid z coordinate

from booms_locations import y_coordinates, z_coordinates, centroid

#Step 2: calculating the moment of inertia of each part around its own centroid:
#and adding the steiner term wrt to the centroid of the crosssection

#Step2a: stiffeners
from booms_locations import boom_area
Izz_stiffeners = np.sum(boom_area * y_coordinates**2)
Iyy_stiffeners = np.sum(boom_area * (z_coordinates-centroid)**2)
#Step2b: Spar
Izz_spar = ((t_sp * h**3)/12)
Iyy_spar = t_sp * h * ((h/2)-centroid)**2
#Step2c: Inclined panel
from booms_locations import l_skin, theta_skin
Izz_plane = (((2*l_skin)**3)* t_sk * (sin(theta_skin))**2)/12 
Iyy_plane = (((2*l_skin)**3)* t_sk * (cos(theta_skin))**2)/12 + 2*l_skin*t_sk*(Ca-centroid)**2
#Step2d: semicircle
centroid_semicircle = pi*h/4

Izz_semicircle = (h/2)**3 * t_sk * pi/2 
Iyy_semicircle = (h/2)**3 * t_sk * pi/2  + pi * (h/2) * t_sk * ((centroid-centroid_semicircle)**2-(h/2-centroid_semicircle)**2)

#Step 3: Calculating finl moment of iertia of the whole cross section

Izz_Aileron = Izz_stiffeners + Izz_spar + Izz_plane + Izz_semicircle
Iyy_Aileron = Iyy_stiffeners + Iyy_spar + Iyy_plane + Iyy_semicircle

print(Izz_Aileron,Iyy_Aileron)
