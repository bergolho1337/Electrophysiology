# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_solution (filename):
    data = np.genfromtxt(filename, delimiter=' ')
    return data

def show_solution (data):
    plt.plot(data[:,0],data[:,1],label="V",c="black")

def write_solution (direction):
    #plt.grid()
    plt.xlabel(u"t (ms)",fontsize=15)
    plt.ylabel(u"V (mV)",fontsize=15)
    #plt.xlim([0,5000])
    #plt.xlim([10000,20000])
    #plt.xlim([20000,30000])
    #plt.xlim([30000,40000])
    #plt.ylim([-90,40])
    #plt.xlim([2400,3000])
    plt.title(u"Action Potential",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    if (direction == "window"):
    	plt.show()
    elif (direction == "pdf"):
	plt.savefig("outputs/output.pdf")
    else:
	print("[-] ERROR! Invalid output option '%s'" % direction)

def main():

    if (len(sys.argv) != 3):
        print("==================================================================")
        print("Usage:> python plot_ap.py <solution_file> <output_direction>")
        print("==================================================================")
	print("<solution_file> = File with the transmembrane potential from the simulation")
	print("<output_direction>:	'window' (Display the figure)")
	print("                         'pdf'	 (Save as a file 'output.pdf')")
	print("==================================================================")
        sys.exit(1)

    solution_filename = sys.argv[1]
    output_direction = sys.argv[2]

    sv = read_solution(solution_filename)

    show_solution(sv)
    write_solution(output_direction)


if __name__ == "__main__":
    main()
