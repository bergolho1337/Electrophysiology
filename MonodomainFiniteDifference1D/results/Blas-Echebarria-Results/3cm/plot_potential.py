import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_potential (sv,output_file):
    plt.grid()
    plt.plot(sv[:,0],sv[:,1],label="Vm",c="black",linewidth=3.0)
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    #plt.xlim(150,560)
    plt.title("Action potential",fontsize=14)
    plt.legend(loc=2,fontsize=14)
    #plt.show()
    plt.savefig(output_file)

def main():

    if (len(sys.argv) != 3):
        print("==========================================================")
        print("Usage:> %s <input_file> <output_file>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_file = sys.argv[1]
	output_file = sys.argv[2]

        data = np.genfromtxt(input_file)

        plot_potential(data,output_file)        
    

if __name__ == "__main__":
    main()
