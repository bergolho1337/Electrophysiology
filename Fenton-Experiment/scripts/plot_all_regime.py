
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def process_data ():
	lmax_arr = [2,3,4,5]
	marks = ['o','D','s','X']
	colors = ['g','b','purple','r']
	labels = ['Normal','Discordant','Corcordant','Block']	
	
	# Plot dummy points for the legend setup
	for i in range(len(marks)):
		plt.scatter(-1,-1,marker=marks[i],color=colors[i],label=labels[i])

	# Plot the data
	for i in range(len(lmax_arr)):
		lmax = i + 2.0
		filename = "Regime/regime-%dcm.dat" % (lmax)
		data = np.genfromtxt(filename)

		for i in range(len(data)):
			bcl, regime = int(data[i][0]), int(data[i][1])
			m = marks[regime]
			c = colors[regime] 		
			plt.scatter(bcl,"%lf" % lmax,marker=m,color=c)
	plt.grid()
	plt.legend(loc=0)
	#plt.show()
	plt.savefig("Regime/regime-all-cables.pdf")

def main():

    if (len(sys.argv) != 1):
        print("==========================================================")
        print("Usage:> %s " % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
 
	process_data()

if __name__ == "__main__":
    main()
