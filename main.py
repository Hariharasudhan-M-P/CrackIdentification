import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
import math as math
plt.style.use('bmh')

#Defining Constants
E  = 200000000000.
m = 7850*.5*.5
pi4 = math.pi*math.pi*math.pi*math.pi
W = .5
L = 3.
b = 1.
I = (.5**4)/12
sinpibl = math.sin(math.pi*b/L)
sin2pibl = math.sin(2*math.pi*b/L)
sin3pibl = math.sin(3*math.pi*b/L)

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
    #return 8*pi4*E*I*.5/m/l4 * (1 + f0*2*W/L*sin2pibl*sin2pibl)
    #return 81 * pi4 * E * I * .5 / m / l4 * (1 + f0 * 2 * W / L * sin3pibl * sin3pibl)
#Mode shape before crack
def phi1sq(x, f0):
    ph = phifactor*(np.sin(np.pi*x/L) - x*b*W/l3*f0*np.pi*np.pi*sinpibl*(1-L/b))
    return ph*ph

#Mode shape after crack
def phi2sq(x, f0):
    ph =  phifactor*(np.sin(np.pi*x/L) + b*W/l2*f0*np.pi*np.pi*sinpibl*(1-x/L))
    return ph*ph

#Mode shape before crack
def phi1sq_2(x, f0):
    ph = phifactor*(math.sin(2*math.pi*x/L) - 4*x*b*W/l3*f0*math.pi*math.pi*sin2pibl*(1-L/b))
    return ph*ph

#Mode shape after crack
def phi2sq_2(x, f0):
    ph =  phifactor*(math.sin(2*math.pi*x/L) + 4*b*W/l2*f0*math.pi*math.pi*sin2pibl*(1-x/L))
    return ph*ph

def phi1sq_3(x, f0):
    ph = phifactor*(math.sin(3*math.pi*x/L) - 9*x*b*W/l3*f0*math.pi*math.pi*sin3pibl*(1-L/b))
    return ph*ph

#Mode shape after crack
def phi2sq_3(x, f0):
    ph =  phifactor*(math.sin(3*math.pi*x/L) + 9*b*W/l2*f0*math.pi*math.pi*sin3pibl*(1-x/L))
    return ph*ph

crack_factor = f(e)
max_kinetic_energy = np.zeros(len(crack_factor)) #kinetic energy
max_potential = u_max(crack_factor)
omega_n = np.zeros(len(crack_factor))
for i in range(len(crack_factor)):
    f0 = crack_factor[i]
    max_kinetic_energy[i] = .5 * m * (si.quad(phi1sq, 0, b, args=(f0))[0] + si.quad(phi2sq, b, L, args=(f0))[0])
    omega_n[i] = math.sqrt(max_potential[i] / max_kinetic_energy[i])

fig,ax1 = plt.subplots(1,1, figsize=(16, 9))
ax1.plot(e, omega_n)
ax1.set_title("Variation of fundamental frequency with crack ratio")
ax1.set_xlabel("Crack Ratio (e)")
ax1.set_ylabel("Natural Frequency $\omega_n$")
#
# fig2 , ax2 = plt.subplots(1,1, figsize=(16, 9))
# ax2.plot(e , max_potential,'-o',label="max potential")
# ax2.plot(e,max_kinetic_energy,'-v',label = "max kinetic")
# ax2.set_title("Variation of Maximum Potential and Kinetic Energy vs crack ratio")
# ax2.set_xlabel("Crack Ratio (e)")
# ax2.set_ylabel("Energies")
# ax2.legend()
#
#
#
# #Mode shape before crack
# def phi1(x, f0):
#     ph = phifactor*(np.sin(np.pi*x/L) - x*b*W/l3*f0*np.pi*np.pi*sinpibl*(1-L/b))
#     return ph
#
# #Mode shape after crack
# def phi2(x, f0):
#     ph =  phifactor*(np.sin(np.pi*x/L) + b*W/l2*f0*np.pi*np.pi*sinpibl*(1-x/L))
#     return ph
#
#
def dphi1(x,f0):
    return  phifactor*(np.pi/L*np.cos(np.pi*x/L) - b*np.pi*np.pi*f0/l3*np.sin(b*np.pi/L)*(1-L/b))

def dphi2(x,f0):
    return  phifactor*(np.pi/L*np.cos(np.pi*x/L) - b*np.pi*np.pi*f0/l3*np.sin(b*np.pi/L))
#
# fig3 = plt.figure(figsize=(16,9))
#
# gs = fig3.add_gridspec(1,4)
# ax = fig3.add_subplot(gs[0, 0:3])
# ay = fig3.add_subplot(gs[0, 3])
#
# ax.set_title("Mode Shape slope vs X for crack ratio = .2")
# ay.set_title("Magnified at 4")
#
# ax.set_xlabel("x")
# ax.set_ylabel("slope")
# ay.set_xlabel("x")
# ay.set_ylabel("slope")
# xx1 = np.linspace(0,4,10)
# xx2 = np.linspace(4,30,70)
# ax.plot(xx1,dphi1(xx1,f(.2)),label="$\\frac{\partial\phi_1(x)}{\partial{x}}$")
# ax.plot(xx2,dphi2(xx2,f(.2)),label="$\\frac{\partial\phi_2(x)}{\partial{x}}$")
# ax.legend(fontsize="x-large")
# ay.plot(xx1,dphi1(xx1,f(.2)),label="$\\frac{\partial\phi_1(x)}{\partial{x}}$")
# ay.plot(xx2,dphi2(xx2,f(.2)),label="$\\frac{\partial\phi_2(x)}{\partial{x}}$")
# ay.legend(fontsize="x-large")
# ay.set_xlim(3,5)
#
#
#
fig4 = plt.figure(figsize=(16,9))

gs = fig4.add_gridspec(1,4)
ax = fig4.add_subplot(gs[0, 0:3])
ay = fig4.add_subplot(gs[0, 3])

ax.set_title("Slope of Mode Shape vs X for crack ratio = .2")
ay.set_title("Magnified at 1")

ax.set_xlabel("x")
ax.set_ylabel("$\phi(x)$")
ay.set_xlabel("x")
ay.set_ylabel("$\phi(x)$")
xx1 = np.linspace(0,1,10)
xx2 = np.linspace(1,3,70)
ax.plot(xx1,dphi1(xx1,f(.2)),label="$\phi_1(x)$")
ax.plot(xx2,dphi2(xx2,f(.2)),label="$\phi_2(x)$")
ax.legend(fontsize="x-large")
ay.plot(xx1,dphi1(xx1,f(.2)),label="$\phi_1(x)$")
ay.plot(xx2,dphi2(xx2,f(.2)),label="$\phi_2(x)$")
ay.legend(fontsize="x-large")
ay.set_xlim(.5,1.5)
#

plt.show()