import mn
import numpy as np
import matplotlib.pyplot as plt
import sys

def simSym(size,steps,sigma,decay):
    a = mn.Forrest(size,decay)
    a.growNetworkSym(sigma)
    for i in range(steps):
        a.collect()
        a.share()
    return a.treeresource

def simR(reps,size,steps,sigma,decay):
    dataS = np.zeros(reps)
    dataP = np.zeros(reps)
    for r in range(reps):
        th = simSym(size,steps,sigma,decay)
        dataS[r] = np.mean(th)
        dataP[r] = np.count_nonzero(th>0.1)/size
    return dataS, dataP

# python sim.py 6 1000 100 3 10

sqrtsize = int(sys.argv[1])   #6
steps = int(sys.argv[2])      #1000
reps = int(sys.argv[3])       #100
decaysteps = int(sys.argv[4])  #3
sigmasteps = int(sys.argv[5]) #10

size = sqrtsize**2
decayrange = np.linspace(0.0,0.002,decaysteps)
sigmarange = np.linspace(0.1,1.0,sigmasteps)

dataS = np.zeros((sigmasteps,decaysteps))
dataP = np.zeros((sigmasteps,decaysteps))

t = 0
for decay in decayrange:
    print("Decay:",decay)
    s = 0
    for sigma in sigmarange:
        print("\tSigma:",sigma)
        dS,dP = simR(reps,size,steps,sigma,decay)
        dataS[s][t] = np.mean(dS)
        dataP[s][t] = np.mean(dP)
        s+=1
    t+=1

np.save("sigmarange.npy",sigmarange)
np.save("resourcedataS.npy",dataS)
np.save("resourcedataP.npy",dataP)

plt.plot(sigmarange,dataS,'o-')
plt.xlabel("Connectivity (sigma)")
plt.ylabel("Average final resources")
plt.show()

plt.plot(sigmarange,dataP,'o-')
plt.xlabel("Connectivity (sigma)")
plt.ylabel("Average trees alive")
plt.show()
