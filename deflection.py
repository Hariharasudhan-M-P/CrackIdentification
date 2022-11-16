import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import math as math
plt.style.use('bmh')

m = 150000
L = 30
E = 50000000000

def eta_uncracked(t):
    return 1.285 * np.sin(.25*np.pi*t)

def eta_cracked(t):
    return -4.1696*np.sin(.25*np.pi*t)

def phi_x(x):
    return np.sqrt(2/(m*L))*np.sin(np.pi*x/L)

#t = [0.33241,1.234235,2.3244,3.434]

t = np.linspace(0,30,50)
fig = plt.figure()

gs = fig.add_gridspec(2,2)
az = fig.add_subplot(gs[1, 0])
ay = fig.add_subplot(gs[1, 1])
ax = fig.add_subplot(gs[0, :])


ax.plot(t,phi_x(4)*eta_uncracked(t),'-o',label="un crcacked x = {x}".format(x = 4))
ax.plot(t,phi_x(4)*eta_cracked(t),'-D',label="cracked x = {x}".format(x = 4))
az.plot(t,phi_x(1)*eta_cracked(t),'-v',label="cracked x = {x}".format(x = 1))
ay.plot(t,phi_x(10)*eta_uncracked(t),'-o',label="uncracked x = {x}".format(x = 10))
ay.plot(t,phi_x(10)*eta_cracked(t),'-D',label="cracked x = {x}".format(x = 10))
az.plot(t,phi_x(10)*eta_uncracked(t),'-v',label="uncracked x = {x}".format(x = 1))
ax.legend()
ay.legend()
az.legend()
ax.set_title("Response for x = 4")
ay.set_title("Response for x > 4")
az.set_title("Response for x < 4")

plt.show()
