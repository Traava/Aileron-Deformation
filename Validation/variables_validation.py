import numpy as np 
#E = 300*10**9      #[Pa]
P = 97.4*10**3     #[N]
theta_max = 28*(np.pi/180) #[rad]
d1 = 1.154/100    #[m]
d3 = 1.840/100    #[m]
n_st = 15         #[-]
w_st = 1.9/100   #[m]
h_st = 1.6/100   #[m]
t_st = 1.2/1000   #[m]
t_sp = 2.8/1000   #[m]
t_sk = 1.1/1000   #[m]
h = 20.5/100      #[m]
xa = 35.0/100     #[m]
x1 = 0.172      #[m]
x2 = 1.211       #[m]
x3 = 2.591        #[m]
la = 2.661        #[m]  
Ca = 0.605        #[m] 
xaI = x2-(xa/2)   #[m]
xaII = x2+(xa/2)  #[m]
#Nx = 41
#Nz = 81
#zhat = 0.0837      #[m] This is just a placeholder guess for now. Real shear center loc needs to be found

G = 27.1e9
E = 72.9e9