import mn
import numpy as np
import matplotlib.pyplot as plt
import sys

g=1.618033
fs = 3.8

sigmarange = np.load("sigmarange.npy")
dataS = np.load("resourcedataS.npy")
dataP = np.load("resourcedataP.npy")

plt.figure(figsize=(fs*g,fs))
plt.plot(sigmarange,dataS,'o-')
plt.title("Role of Nutrient Decay")
plt.xlabel("Connectivity (sigma)")
plt.ylabel("Average final resources")
plt.legend(["0.000","0.001","0.002"], frameon=False)
plt.savefig("avgfinalres.pdf")
plt.show()

plt.figure(figsize=(fs*g,fs))
plt.plot(sigmarange,dataP,'o-')
plt.title("Role of Nutrient Decay")
plt.xlabel("Connectivity (sigma)")
plt.ylabel("Average trees alive")
plt.legend(["0.000","0.001","0.002"], frameon=False)
plt.savefig("avgtreealive.pdf")
plt.show()
