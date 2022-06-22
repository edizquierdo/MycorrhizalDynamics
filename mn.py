import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

class Forrest:

    def __init__(self, size, decay):
        self.treenumber = size
        # density = 0 is a probability that any two trees will be connected by a fungus
        self.adjmat = np.zeros((size,size))
        #print("MN adjacency matrix:\n",self.adjmat)
        d = int(np.sqrt(size))
        self.treeposx = np.array([np.linspace(-1,1,d)]*d).flatten()
        self.treeposy = np.array([np.linspace(-1,1,d)]*d).T.flatten()
        #print("Physical position of trees:\n",self.treeposx,self.treeposy)
        self.resourceposx = 0.0 #np.random.random()
        self.resourceposy = 0.0 #np.random.random()
        #print("Position of resource:",self.resourcepos)
        self.treeresource = np.zeros(size) #np.array([0.]*size) #np.random.random(size)
        #print("Starting resources of trees:",self.treeresource)
        self.sigma = 0.2     # the broadness of the normal XXX Change that
        self.mu = 0.0        # the location o the normal distribution of the nutrient
        self.threshold = 1.0
        self.gain = 0.01
        self.decay = decay  #0.001
        self.transportcost = 0.1
        self.maxreplimit = size

    def growNetworkSym(self,sigma):
        distances = np.zeros((self.treenumber,self.treenumber))
        for i in range(self.treenumber):
            for j in range(i+1,self.treenumber):
                # Calculate distance between tree i and tree j
                dist = np.sqrt((self.treeposx[i]-self.treeposx[j])**2+(self.treeposy[i]-self.treeposy[j])**2)
                # Collect proportional to distance (assume some Gaussian or Normal Distribution)
                prob = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(dist**2) / (2*(sigma**2)))
                # Determine if fungal connection is made or not based on proximity
                if np.random.random() < prob:
                    self.adjmat[i,j] = 1
                    self.adjmat[j,i] = 1

    def growNetworkAsym(self,sigma):
        distances = np.zeros((self.treenumber,self.treenumber))
        for i in range(self.treenumber):
            for j in range(self.treenumber):
                # Calculate distance between tree i and tree j
                dist = np.sqrt((self.treeposx[i]-self.treeposx[j])**2+(self.treeposy[i]-self.treeposy[j])**2)
                # Collect proportional to distance (assume some Gaussian or Normal Distribution)
                prob = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(dist**2) / (2*(sigma**2)))
                # Determine if fungal connection is made or not based on proximity
                if np.random.random() < prob:
                    self.adjmat[i,j] = 1

    def showTreePos(self):
        plt.plot(self.treeposx,self.treeposy,'go')
        for i in range(self.treenumber):
            for j in range(self.treenumber):
                if self.adjmat[i,j] == 1:
                    plt.plot([self.treeposx[i],self.treeposx[j]],[self.treeposy[i],self.treeposy[j]],'-')
        plt.show()

    def neighbors(self,i):
        return [j for j in range(self.treenumber) if self.adjmat[i, j]]

    def collect(self):
        # Calculate distance between tree i and resource
        dist = np.sqrt((self.treeposx-self.resourceposx)**2+(self.treeposy-self.resourceposy)**2)
        # Collect proportional to distance (assume some Gaussian or Normal Distribution)
        nutrients = (1/(self.sigma*np.sqrt(2*np.pi)))*np.exp(-((dist - self.mu)**2) / (2*(self.sigma**2)))
        # Add it to current tree
        self.treeresource += self.gain * nutrients - self.decay
        # Clip to zero; no negative nutrients -- also, if a tree hits zero, should it disappear?
        self.treeresource = np.clip(self.treeresource,0.0,1.5)

    def share(self):
        reps = 0
        # 1. Repeat sharing until there's no tree with more nutrient than the treshold
        while (np.max(self.treeresource) > self.threshold) and (reps < self.maxreplimit):
            # 2. Find the tree that's highest above threshold
            i = np.argmax(self.treeresource)
            # 3. Calculate its excess
            excess = self.treeresource[i] - self.threshold
            # 4. Remove its excess
            self.treeresource[i] = self.treeresource[i] - excess
            # 5. Find its neighbors
            neighbors = self.neighbors(i)
            if (len(neighbors)>0):
                # 6. Calculate share evenly (there's % transport cost)
                excessshare = (excess/len(neighbors))*(1-self.transportcost)
                # 7. Share extra resources to neighbors
                for j in neighbors:
                    self.treeresource[j] += excessshare
            # Keep track of how many times you've done this
            reps += 1
