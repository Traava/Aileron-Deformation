from ShearDistribution import *
import numpy as np
from MatrixSystem import T, Vy, Vz
import matplotlib.pyplot as plt

Shear_Distribution = np.array([])


def ShearPlots(n):
    
    x = np.linspace(0,la,n) 
    
    max_q = np.array([0,0])
    
    for i in range(n):     
        
        span_location = i*(la/n)
        
        Sy = Vy(x)[i]
        Sz = Vz(x)[i]
        Ts  = T(x)[i]
        
        qb1, qb2, qb3, qb4, qb5, qb6 = ShearDistribution(Ts, Sy, Sz)
        tau1 = qb1/t_sk
        tau2 = qb2/t_sp
        tau3 = qb3/t_sk
        tau4 = qb3/t_sk
        tau5 = qb5/t_sp
        tau6 = qb6/t_sk
        
        tot_shear = np.hstack((tau1,tau2,tau3,tau4,tau5,tau6))
        
        max_q_x = np.array([span_location, np.max(tot_shear)])
        
        max_q = np.vstack((max_q, max_q_x))
        
    #print(np.max(np.amax(max_q, axis = 1)))
    Max_location = np.where(max_q == np.max(np.amax(max_q, axis = 0)))
    Max_loc = np.asscalar(Max_location[0])
    loc = (Max_loc-1)*la/n
    

    print((np.max(np.amax(max_q, axis = 1)))*(10**-6))
    print(loc)
    #print(Max_loc)
    #print(np.shape(tot_shear))
    
    
    """
    
    Syy = Vy(x)[Max_loc]
    Szz = Vz(x)[Max_loc]
    Tss  = T(x)[Max_loc]
        
    qb11, qb22, qb33, qb44, qb55, qb66 = ShearDistribution(Tss, Syy, Szz)
    tau11 = qb11/t_sk
    tau22 = qb22/t_sp
    tau33 = qb33/t_sk
    tau44 = qb44/t_sk
    tau55 = qb55/t_sp
    tau66 = qb66/t_sk
     
    tot_shearr = np.hstack((tau11,tau22,tau33,tau44,tau55,tau66))
    
    #z = np.linspace(0, l_sk*2 + pi*h/2, 700)
    #plt.plot(z, np.hstack((tau1,tau3,tau4,tau6)))
    
    z_coords = np.zeros(800)
    y_coords = np.zeros(800)
    #region1.1
    for i in range(50):
        sep1 = theta_separation/50
        y_coords[i] = -r*sin(sep1*i)
        z_coords[i] = -centroid + r*(1-cos(sep1*i))
    #region1.2  
    for i in range(50):
        sep2 = (pi*r/2 - theta_separation)/50
        y_coords[i + 50] = -r*sin(sep2*i + theta_separation)
        z_coords[i + 50] = -centroid + r*(1-cos(sep2*i + theta_separation))
     #region2
    for i in range(50):
        sep3 = r/50
        z_coords[i + 100] = r - centroid
        y_coords[i + 100] = -sep3*i
    #region3.1
    for i in range(50):
        sep4 = (booms_separation*2 - pi*r/2)/50
        y_coords[i + 150] = -h/2 + sep4*i*sin(theta_skin)
        z_coords[i + 150] = -h/2 + sep4*i*cos(theta_skin)
    #region3.2 
    for i in range(150):
        sep5 = booms_separation/50
        y_coords[i + 200] = -h/2 + (2*booms_separation -  pi*r/2 + sep5*i)*sin(theta_skin) 
        z_coords[i + 200] = -h/2 + (2*booms_separation -  pi*r/2 + sep5*i)*cos(theta_skin) 
    for i in range (50):
        sep6 = booms_separation/100
        y_coords[i + 350] = -h/2 + (l_skin - 0.5*booms_separation + sep6*i)*sin(theta_skin)
        z_coords[i + 350] = ((l_skin - 0.5*booms_separation + sep6*i)*cos(theta_skin) + h/2) - centroid
    #for i in range(400):
     #   y_coords[i+400] = -y_coords[400-i]
      #  z_coords[i+400] = z_coords[400-i]
        
    
    for i in range (50):
        sep7 = booms_separation/100
        y_coords[i+400] = sep7*i*sin(theta_skin)
        z_coords[i+400] = (Ca-centroid) - sep7*i*cos(theta_skin)
    for i in range(150):
        sep8 = booms_separation/50
        y_coords[i+450] = (booms_separation/2 + sep8*i)*sin(theta_skin)
        z_coords[i+450] = (Ca-centroid) - (booms_separation/2 + sep7*i)*cos(theta_skin)
    for i in range(50):
        sep9 = sep4
        y_coords[i+600] = (booms_separation*3.5 + sep9*i)*sin(theta_skin)
        z_coords[i+600] = (Ca-centroid) - (booms_separation*3.5 + sep7*i)*cos(theta_skin)
    for i in range(50):
        sep10 = sep3
        z_coords[i + 650] = r - centroid
        y_coords[i + 650] = sep9*i
    for i in range(50):
        sep11 = sep2
        y_coords[i + 700] = r*sin(pi/2 - sep11*i)
        z_coords[i + 700] = -centroid + r*cos(pi/2 -sep11*i)
    for i in range(50):
        sep12 = sep1
        y_coords[i+750] = r*sin(theta_separation - sep11*i)
        z_coords[i+750] = -centroid + r*cos(theta_separation -sep11*i)
       
    
    
    
    fig1 = plt.figure()
    ax1 = plt.axes()
    #ax1.invert_yaxis()
    img1 = plt.scatter(z_coords, y_coords, c = tot_shear, cmap ='RdYlGn_r')
    #img2 = plt.scatter(z_coord_sp, y_coord_sp, c = Stress_sp[20], cmap ='RdYlGn_r') 
    
    cbar1 = fig1.colorbar(img1)
    cbar1.set_label(r'[$\frac{N}{m^2}$]')
    ax1.set_title("Bending Stress Distribution")
    ax1.set_xlabel('z [m]')
    ax1.set_ylabel('y [m]')
        
    #print(z_coords)
    #print(y_coords)

        
    
    
    
    
    
    
    
    #return Max_location
    #print(np.shape(Max_location))
    #loc = np.asscalar(Max_location[0][0])
    
    #print('Max Shear Stress occurs at x= ', Max_location*la/n )
    
    #return loc*la/n
    """
    
        
        
        
        
        
    