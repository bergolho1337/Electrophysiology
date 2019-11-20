import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from termcolor import colored

def print_header (title,colorname):
	print(colored(title,colorname))
	print("     "),
	for bcl in range(280,145,-5):
		print("%d " % bcl),
	print(" {ms}")

def print_result (filenames,title,colorname):
	
	print_header(title,colorname)

	cable_length = 5.0

	for filename in filenames:
		result = []

		with open(filename) as fp:
			line = fp.readline()
			while line:
				tokens = line.split()
				state = tokens[1]
			
				if (state == "No"):
					result.append("N")
				elif(state == "Discordant"):
					result.append("D")
				elif(state == "Concordant"):
					result.append("C")
				elif(state == "Block"):
					result.append("B")

				line = fp.readline()

			print("%.1lfcm  " % cable_length),
			for i in range(len(result)):
				state = result[i]
				if (state == "N"):
					print(colored("%s   " % result[i], "green")),
				elif(state == "D"):
					print(colored("%s   " % result[i], "yellow")),
				elif(state == "C"):
					print(colored("%s   " % result[i], "blue")),
				elif(state == "B"):
					print(colored("%s   " % result[i], "red")),
			print("")
		
		cable_length = cable_length - 0.5

	print(colored("N -> No alternans","green"))
	print(colored("D -> Discordant alternans","yellow"))
	print(colored("C -> Concordant alternans","blue"))
	print(colored("B -> Block propagation","red"))


def main():

	if (len(sys.argv) != 1):
		print("==========================================================")
		print("Usage:> %s" % (sys.argv[0]))
		print("==========================================================")
		sys.exit(1)

	print("******************************************************************************************************************************************************")
	karma_filenames = ["karma-results/alternans-5cm.dat","karma-results/alternans-4.5cm.dat","karma-results/alternans-4cm.dat","karma-results/alternans-3.5cm.dat","karma-results/alternans-3cm.dat","karma-results/alternans-2.5cm.dat","karma-results/alternans-2cm.dat","karma-results/alternans-1.5cm.dat","karma-results/alternans-1cm.dat"]
	karma_title = "[Result] Spatiotemporal control of cardiac alternans (Blas Echebarria and Alain Karma) 2002" 
	print_result(karma_filenames,karma_title,"magenta")
	print("******************************************************************************************************************************************************")

	print("******************************************************************************************************************************************************")
	alejandro_filenames = ["alejandro-results/alternans-5cm.dat","alejandro-results/alternans-4.5cm.dat","alejandro-results/alternans-4cm.dat","alejandro-results/alternans-3.5cm.dat","alejandro-results/alternans-3cm.dat","alejandro-results/alternans-2.5cm.dat","alejandro-results/alternans-2cm.dat","alejandro-results/alternans-1.5cm.dat","alejandro-results/alternans-1cm.dat"]
	alejandro_title = "[Result] Model-based control of cardiac alternans in Purkinje fibers (Alejandro Garzon and Roman O. Grigoriev and Flavio Fenton) 2011"
	print_result(alejandro_filenames,alejandro_title,"cyan")
	print("******************************************************************************************************************************************************")

	print("******************************************************************************************************************************************************")
	lucas_filenames = ["lucas-results/alternans-5cm.dat","lucas-results/alternans-4.5cm.dat","lucas-results/alternans-4cm.dat","lucas-results/alternans-3.5cm.dat","lucas-results/alternans-3cm.dat","lucas-results/alternans-2.5cm.dat","lucas-results/alternans-2cm.dat","lucas-results/alternans-1.5cm.dat","lucas-results/alternans-1cm.dat"]
	lucas_title = "[Result] Lucas Berg (Homogenous model) 2019"
	print_result(lucas_filenames,lucas_title,"yellow")
	print("******************************************************************************************************************************************************")

if __name__ == "__main__":
    main()
