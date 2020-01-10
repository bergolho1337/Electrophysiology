#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import numpy as np

def main():

	if len(sys.argv) != 1:
		print("-------------------------------------------------------------------------")
                print("Usage:> python %s " % sys.argv[0])
                print("-------------------------------------------------------------------------")
                return 1

    	data = np.genfromtxt("tables/output_table.dat")

	# Sorting by the 3 column, which is the ATPI value
	sorted_data = np.sort(data.view("f8,f8,f8,f8,f8,f8,f8,f8,f8,f8"),order=["f3"],axis=0).view(np.float)

	output_file = open("tables/sorted_table.dat","w")
	for i in range(len(sorted_data)):
		output_file.write("%g %g %g %g %g %g %g %g %g %g\n" % (sorted_data[i][0],sorted_data[i][1],sorted_data[i][2],sorted_data[i][3],sorted_data[i][4],sorted_data[i][5],sorted_data[i][6],sorted_data[i][7],sorted_data[i][8],sorted_data[i][9]))
	output_file.close()

if __name__ == "__main__":
    	main()
