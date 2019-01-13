import os
import sys
import numpy as np
from matplotlib import pyplot

def main ():
    data = np.genfromtxt(open("output.dat","r"))
    
    pyplot.clf()
    pyplot.title("Action Potential")
    pyplot.xlabel("Time (s)")
    pyplot.ylabel("Potential (mV)")
    pyplot.plot(data[:,0],data[:,1],label="V",linewidth=2,color="black")
    pyplot.grid()
    pyplot.legend(loc=0,fontsize=15)
    pyplot.savefig("v.pdf")

if __name__ == "__main__":
    main()
