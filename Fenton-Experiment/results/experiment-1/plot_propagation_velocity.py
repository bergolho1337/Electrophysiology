import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_propagation_velocity (positions,velocities,cable_length,sigma): 
    # Take the first cell out
    positions = positions[1:]
    velocities = velocities[1:]

    # Plot the data
    plt.plot(positions,velocities,label=r"$\sigma$ = %lf" % (sigma),linewidth=3.0)

    #plt.savefig("ap.pdf")

def calculate_cell_positions (cable_length,dx):
    num_cells = int(cable_length/dx)
    offset = num_cells/5
    
    cell_positions = []
    for i in range(5):
        cell_position = i*offset*dx
        cell_positions.append(cell_position)
    
    return cell_positions

def get_cell_velocities_from_file (input_file):
    data = np.genfromtxt(input_file)
    velocities = data[:,1]

    return velocities

def main():

    if (len(sys.argv) != 3):
        print("==========================================================")
        print("Usage:> %s <input_folder> <dx>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_folder = sys.argv[1]
        dx = float(sys.argv[2])

        cable_lengths = np.arange(10000,55000,5000)
        #cable_lengths = np.arange(10000,15000,5000)
        sigmas = np.arange(0.0001,0.001,0.0001)

        for i in range(len(cable_lengths)):
            for j in range(len(sigmas)):
                cable_length = cable_lengths[i]
                sigma = sigmas[j]

                input_file = "%s/cable-%.1lfum-sigma-%.6lf/propagation_velocity.dat" % (input_folder,cable_length,sigma)
                cell_positions = calculate_cell_positions(cable_length,dx)
                cell_velocities = get_cell_velocities_from_file(input_file)

                plot_propagation_velocity(cell_positions,cell_velocities,cable_length,sigma)

            plt.grid()
            plt.xlabel("cell position (um)",fontsize=15)
            plt.ylabel("v (m/s)",fontsize=15)
            plt.xlim(0,cable_lengths[-1])
            #plt.ylim(0,4)
            plt.title("Velocity x Cell position",fontsize=14)
            plt.legend(loc=0,fontsize=14)
            plt.show()    
            plt.clf()

if __name__ == "__main__":
    main()
