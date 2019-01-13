# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_solution (filename):
    data = np.genfromtxt(filename, delimiter=' ')
    return data

def show_solution (data):
    plt.scatter(data[:,0],data[:,1],label="V",c="black",marker='o')

def write_solution ():
    plt.grid()
    plt.xlabel(u"BCL (ms)",fontsize=15)
    plt.ylabel(u"APD (ms)",fontsize=15)
    #plt.xlim([0,5000])
    #plt.xlim([10000,20000])
    #plt.xlim([20000,30000])
    #plt.xlim([30000,40000])
    plt.title(u"BCL x APD",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()
    #plt.savefig("output.pdf")

def main():

    if (len(sys.argv) != 2):
        print("========================================================")
        print("Usage:> python bcl_apd_plot.py <bcl_apd_file>")
        print("========================================================")
        sys.exit(1)

    solution_filename = sys.argv[1]

    sv = read_solution(solution_filename)
    
    show_solution(sv)
    write_solution()


if __name__ == "__main__":
    main()