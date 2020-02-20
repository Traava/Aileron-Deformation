import numpy as np 
E = 300*10**9      #[Pa]
P = 49.2*10**3     #[N]
theta_max = 30*(np.pi/180) #[rad]
d1 = 0.389/100    #[m]
d3 = 1.245/100    #[m]
n_st = 11         #[-]
w_st = 1.7/1000   #[m]
h_st = 1.3/1000   #[m]
t_st = 1.2/1000   #[m]
t_sp = 2.4/1000   #[m]
t_sk = 1.1/1000   #[m]
h = 16.1/100      #[m]
xa = 24.5/100     #[m]
x1 = 0.125        #[m]
x2 = 0.498        #[m]
x3 = 1.494        #[m]
la = 1.611        #[m]  #1.611
Ca = 0.505        #[m]  #0.505
xaI = x2-(xa/2)   #[m]
xaII = x2+(xa/2)  #[m]
Nx = 41
Nz = 81
zhat = Ca/2       #[m] This is just a placeholder guess for now. Real shear center loc needs to be found


G = 27.1e9
E = 72.9e9