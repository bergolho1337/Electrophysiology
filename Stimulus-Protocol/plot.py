import os
import sys
import numpy as np
from matplotlib import pyplot

def plotSolution (tmin,tmax):
    nameStates = ["V","h"]
    data = np.genfromtxt(open("solution.dat","r"))
    for i in range(len(nameStates)):
        pyplot.clf()
        pyplot.title(nameStates[i] + " x Time")
        pyplot.xlabel("Time")
        pyplot.ylabel(nameStates[i])
        pyplot.xlim(tmin,tmax)
        pyplot.plot(data[:,0],data[:,i+1],label=nameStates[i],linewidth=2,color="black")
        pyplot.grid()
        pyplot.legend(loc=0,fontsize=15)
        #pyplot.savefig(nameStates[i] + ".pdf")
	    #pyplot.savefig("Output/" + nameStates[i] + ".png")
        pyplot.show()

def main ():
    plotSolution(4000,8000)

if __name__ == "__main__":
    main()
