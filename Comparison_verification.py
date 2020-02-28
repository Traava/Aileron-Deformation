from MatrixSystem import Sz_num, Sy_num, Mz_num, My_num, T_num, v_num, w_num, phi_num,xx
from main import aileron
import matplotlib.pyplot as plt
import numpy as np
from variables import *


x = np.linspace(0,la,num = 1000)
Sy_ver = aileron.Sy(x)
Sz_ver = aileron.Sz(x)
My_ver = aileron.My(x)
Mz_ver = aileron.Mz(x)
T_ver = aileron.T(x)
v_ver,w_ver,phi_ver = aileron.eval(x)

shift = phi_num[0]+phi_ver[0]
offset = 0.05
plt.figure()

plt.subplot(221)
plt.plot(xx,v_num, label = "Numerical model",color = 'black')
plt.plot(x,-v_ver, label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("v [m]").set_size(14)
plt.legend()
plt.grid()

plt.subplot(222)
plt.plot(xx,w_num,label = "Numerical model",color = 'black')
plt.plot(x,-w_ver, label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("w [m]").set_size(14)
plt.legend()
plt.grid()

plt.subplot(223)
plt.plot(xx,phi_num-np.ones(len(phi_num))*shift,label = "Numerical model",color = 'black')
plt.plot(x,-phi_ver,label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("$\phi$ [rad]").set_size(14)
plt.legend()
plt.grid()

plt.subplot(224)
plt.plot(xx,T_num/1000,label = "Numerical model",color = 'black')
plt.plot(x,-T_ver/1000,label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("T [kNm]").set_size(14)
plt.legend()
plt.grid()


plt.show()

plt.figure()

plt.subplot(221)
plt.plot(xx,Sy_num/1000,label = "Numerical model",color = 'black')
plt.plot(x,-Sy_ver/1000,label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.ylim(ymin = -30)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("$S_{y}$ [kN]").set_size(14)
plt.legend()
plt.grid()

plt.subplot(222)
plt.plot(xx,Sz_num/1000,label = "Numerical model",color = 'black')
plt.plot(x,-Sz_ver/1000,label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.ylim(ymax = 200)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("$S_{z}$ [kN]").set_size(14)
plt.legend()
plt.grid()

plt.subplot(223)
plt.plot(xx,Mz_num/1000,label = "Numerical model",color = 'black')
plt.plot(x,-Mz_ver/1000,label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("$M_{z}$ [kNm]").set_size(14)
plt.legend()
plt.grid()

plt.subplot(224)
plt.plot(xx,My_num/1000,label = "Numerical model",color = 'black')
plt.plot(x,-My_ver/1000,label = "Verification model", color = 'red',ls =  '--')
plt.xlim(xmin = 0-offset, xmax = la+offset)
plt.xlabel("x [m]").set_size(14)
plt.ylabel("$M_{y}$ [kNm]").set_size(14)
plt.legend()
plt.grid()


plt.show()


print("Maximum torque ver:", max(-T_ver))
print("Maximum torque num:", max(T_num))

print("Maximum My ver:", min(-My_ver))
print("Maximum My num:", min(My_num))

print("Maximum Mz ver:", max(-Mz_ver))
print("Maximum Mz num:", max(Mz_num))

print("Maximum Sy ver:", max(-Sy_ver[200:]))
print("Maximum Sy num:", max(Sy_num))

print("Maximum Sz ver:", min(-Sz_ver[100:]))
print("Maximum Sz num:", min(Sz_num))



