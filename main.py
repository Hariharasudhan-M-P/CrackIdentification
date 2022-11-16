import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import math as math
plt.style.use('bmh')

#Defining Constants
E  = 50000000000
m = 150000.
pi4 = math.pi*math.pi*math.pi*math.pi
W = 1.
L =30.
b = 4.
I = 0.05
sinpibl = math.sin(math.pi*b/L)
phifactor = math.sqrt(2/(m*L))
l4 = L*L*L*L
l3 = L*L*L
l2 = L*L
#Don't wana use pows in ilteration so making them constants
e = np.linspace(0,0.9,50) #crack ratio

#F
def f(e):
    return 2*(e/(1-e))*(e/(1-e))*(5.93 - 19.69*e + 37.14*e*e - 35.84*e*e*e + 13.12*e*e*e*e)

#Potential Energy
def u_max(f0):
    return pi4*E*I*.5/m/l4 * (1 + f0*2*W/L*sinpibl*sinpibl)

#Mode shape before crack
def phi1(x,f0):
    ph = phifactor*(math.sin(math.pi*x/L) - x*b*W/l3*f0*math.pi*math.pi*sinpibl*(1-L/b))
    return ph*ph

#Mode shape after crack
def phi2(x,f0):
    ph =  phifactor*(math.sin(math.pi*x/L) + b*W/l2*f0*math.pi*math.pi*sinpibl*(1-x/L))
    return ph*ph

crack_factor = f(e)
max_kinetic_energy = np.zeros(len(crack_factor)) #kinetic energy
max_potential = u_max(crack_factor)
omega_n = np.zeros(len(crack_factor))
for i in range(len(crack_factor)):
    f0 = crack_factor[i]
    max_kinetic_energy[i] = .5 * m * (si.quad(phi1, 0, b, args=(f0))[0] + si.quad(phi2, b, L, args=(f0))[0])
    omega_n[i] = math.sqrt(max_potential[i] / max_kinetic_energy[i])

fig,ax1 = plt.subplots(1,1, figsize=(16, 9))
ax1.plot(e, omega_n)
ax1.set_title("Variation of fundamental frequency with crack ratio")
ax1.set_xlabel("Crack Ratio (e)")
ax1.set_ylabel("Natural Frequency $\omega_n$")

fig2 , ax2 = plt.subplots(1,1, figsize=(16, 9))
ax2.plot(e , max_potential,'-o',label="max potential")
ax2.plot(e,max_kinetic_energy,'-v',label = "max kinetic")
ax2.set_title("Variation of Maximum Potential and Kinetic Energy vs crack ratio")
ax2.set_xlabel("Crack Ratio (e)")
ax2.set_ylabel("Energies")
ax2.legend()


def dphi1(x,f0):
    return (1/1500)*(np.pi/30*np.cos(np.pi*x/30) - 4*np.pi*np.pi*f0/27000*np.sin(4*np.pi/30)*(1-30/4))

def dphi2(x,f0):
    return (1/1500)*(np.pi/30*np.cos(np.pi*x/30) - 4*np.pi*np.pi*f0/27000*np.sin(4*np.pi/30))

fig3 = plt.figure(figsize=(16,9))

gs = fig3.add_gridspec(1,4)
ax = fig3.add_subplot(gs[0, 0:3])
ay = fig3.add_subplot(gs[0, 3])

ax.set_title("Mode Shape slope vs X for crack ratio = .2")
ay.set_title("Magnified at 4")

ax.set_xlabel("x")
ax.set_ylabel("slope")
ay.set_xlabel("x")
ay.set_ylabel("slope")
xx1 = np.linspace(0,4,10)
xx2 = np.linspace(4,30,70)
ax.plot(xx1,dphi1(xx1,f(.2)),label="$\\frac{\partial\phi_1(x)}{\partial{x}}$")
ax.plot(xx2,dphi2(xx2,f(.2)),label="$\\frac{\partial\phi_2(x)}{\partial{x}}$")
ax.legend(fontsize="x-large")
ay.plot(xx1,dphi1(xx1,f(.2)),label="$\\frac{\partial\phi_1(x)}{\partial{x}}$")
ay.plot(xx2,dphi2(xx2,f(.2)),label="$\\frac{\partial\phi_2(x)}{\partial{x}}$")
ay.legend(fontsize="x-large")
ay.set_xlim(3,5)

plt.show()


