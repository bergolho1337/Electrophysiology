# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_solution (filename):
    data = np.genfromtxt(filename, delimiter=' ')
    return data

def show_solution (data_1,data_2):
    plt.plot(data_2[:,0],data_2[:,1],label="Full fibrotic cell",c="blue",linewidth=1.0)
    plt.plot(data_1[:,0],data_1[:,1],label="Healthy cell",c="red",linewidth=1.0)

def write_solution (direction):
    #plt.grid()
    plt.xlabel(u"t (ms)",fontsize=15)
    plt.ylabel(u"V (mV)",fontsize=15)
    #plt.xlim([11300,12000])
    plt.xlim([-10,1300])
    #plt.xlim([20000,30000])
    #plt.xlim([30000,40000])
    #plt.ylim([-90,40])
    #plt.xlim([2400,3000])
    plt.title(u"Action Potential - INaCa = 100%",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    if (direction == "window"):
    	plt.show()
    elif (direction == "pdf"):
	plt.savefig("outputs/output.pdf")
    else:
	print("[-] ERROR! Invalid output option '%s'" % direction)

def main():

    if (len(sys.argv) != 4):
        print("==================================================================")
        print("Usage:> python plot_ap.py <solution_file_1> <solution_file_2> <output_direction>")
        print("==================================================================")
	print("<solution_file> = File with the transmembrane potential from the simulation")
	print("<output_direction>:	'window' (Display the figure)")
	print("                         'pdf'	 (Save as a file 'output.pdf')")
	print("==================================================================")
        sys.exit(1)

    solution_filename_1 = sys.argv[1]
    solution_filename_2 = sys.argv[2]
    output_direction = sys.argv[3]

    sv_1 = read_solution(solution_filename_1)
    sv_2 = read_solution(solution_filename_2)

    show_solution(sv_1,sv_2)
    write_solution(output_direction)


if __name__ == "__main__":
    main()
