import matplotlib.pyplot as plt
import math as math
import numpy as np
from tqdm import *
from subprocess import call
import os

def writeMedium(N,length,T,XH2O,XCO2):

    f = open('Medium', 'w')
    for i in range(0,N):
        line = str(length[i]) + " " + str(T[i]) + " " + str(XH2O[i]) + " " + str(XCO2[i]) + "\n"
        f.write(line)

    f.close()

# save directory

nameSave = "mixture-T-X-pw1/"

if not os.path.exists(nameSave):
        os.makedirs(nameSave)

# definition of medium

S = 1
ncol = 50
Dx = S/ncol
x = np.arange(0,S+0.0001,Dx)
N = len(x)
xd =np.array([x[i] + (x[i+1]-x[i])/2 for i in range(0,N-1)])

#T = 1250-750*np.cos(2*math.pi*xd/S)
T = 1500*np.exp(-(xd-S/2.)**2/(2*0.01))+300
#T = 400+1400*np.sin(2*math.pi*xd/S)**2
#T = np.ones(ncol)*1000
Tw = 300

plt.plot(xd,T)
plt.xlabel("X [m]")
plt.ylabel("Temperature [K]")
plt.savefig(nameSave+"T.eps", format="eps", dpi=100)
plt.show()
plt.clf()

#T = np.ones(ncol)*1600

length = np.ones(ncol)*Dx*100

#XH2O = np.ones(ncol)*0.6
#XCO2 = np.ones(ncol)*0.3
XCO2 = np.exp(-(xd-S/2.)**2/(2*0.01))*0.6
XH2O = XCO2
XN2 = 1-XCO2-XH2O

plt.plot(xd, XCO2, label="CO2")
plt.plot(xd, XH2O, label="H2O")
plt.plot(xd, XN2, label="XN2")
plt.legend()
plt.xlabel("X [m]")
plt.ylabel("$\chi$ [K]")
plt.savefig(nameSave+"X.eps", format="eps", dpi=100)
plt.show()
plt.clf()


writeMedium(ncol,length,T,XH2O,XCO2)

# run SNB
#

print("--- running SNB ---")

call(["./snb", str(Tw)])

#Intensity

print("--- reading intensity ---")
data = np.loadtxt("I")


print("--- calculating flux ---")
# DOM
#
Ip = data[:,1]
Im = -np.flip(Ip,axis=0)

#flux

q = Ip + Im

plt.plot(x,q, marker = 'o')
plt.show()
plt.clf()


print("--- calculating radiative source term ---")


# source term

divq = [-(q[i+1]-q[i])/(Dx)/1000 for i in range(0,N-1)]


# moving simulations results to directory

os.rename("I", nameSave+"I")
os.rename("Inu", nameSave+"Inu")
os.rename("Medium", nameSave+"Medium")

# save

divqA = np.array(divq)
xA = np.array(xd)
save = np.c_[xA,divqA]
np.savetxt(nameSave+"DivQ",save)

#figures

dataWSGGC=np.loadtxt("../mixture-T-X/1D-0.1/coelho3D_uniformCase1_bandC/postProcessing/divRad/1e-07/line_dQrad.xy")
dataWSGGS=np.loadtxt("../mixture-T-X/1D-0.1/coelho3D_uniformCase1_bandSmith/postProcessing/divRad/1e-07/line_dQrad.xy")
dataWSGGJ=np.loadtxt("../mixture-T-X/1D-0.1/coelho3D_uniformCase1_bandJ/postProcessing/divRad/1e-07/line_dQrad.xy")
dataGrey=np.loadtxt("../mixture-T-X/1D-0.1/grey/postProcessing/divRad/1e-07/line_dQrad.xy")

plt.plot(dataWSGGJ[:,0],dataWSGGJ[:,1]/1000,label="WSGG - Johansson")
plt.plot(dataWSGGC[:,0],dataWSGGC[:,1]/1000,label="WSGG - Cassol")
plt.plot(dataWSGGS[:,0],dataWSGGS[:,1]/1000,label="WSGG - Smith")
plt.plot(dataGrey[:,0],dataGrey[:,1]/1000,label="Grey Model")

plt.plot(xd,divqA,'k', marker='o', label='SNB')

plt.xlabel("X [m]")
plt.ylabel("-divQ$_{rad}$ [kW/m$^{-3}$]")

plt.xlim(0,S)
#plt.ylim(-550,-200)

plt.grid()
plt.legend()

plt.savefig(nameSave+"divQ.eps", format="eps", dpi=100)

plt.show()

