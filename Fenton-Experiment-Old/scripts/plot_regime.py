# To use the 'plot_regime.py' first run the 'calc_apd_bcl_diagram.py' script and 
# copy the last cell output to the Regime folder


import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def process_data (data):
	regime = []
	for i in range(len(data)):
		if (data[i][1] == 0.0 or data[i][2] == 0.0):	
			regime.append(0)
		elif (data[i][1] - data[i][2] == 0.0):
			regime.append(1)
		else:
			regime.append(2)
	return regime

def plot_regime (regs,bcls,lmax):
	marks = ['X','o','D']
	colors = ['r','g','b']	
	for i in range(len(bcls)):
		m = marks[regs[i]]
		c = colors[regs[i]] 		
		plt.scatter(bcls[i],"%lf" % lmax,marker=m,color=c)
	plt.grid()
	#plt.savefig("Regime/regime_over_cable.pdf")	
	plt.show()

def main():

    if (len(sys.argv) != 3):
        print("==========================================================")
        print("Usage:> %s <input_file> <lmax>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_file = sys.argv[1]
	lmax = float(sys.argv[2])

	data = np.genfromtxt(input_file) 

	regime = process_data(data)

	plot_regime(regime,data[:,0],lmax)
    

if __name__ == "__main__":
    main()
