import sys
import subprocess
import time
import numpy as np

def main ():
    if ( len(sys.argv) != 4 ):
        print("-------------------------------------------------------------------------------------------------------------")
        print("Usage:> python %s <input_filename_1> <input_filename_2> <BCL>")
        print("-------------------------------------------------------------------------------------------------------------")
        print("<input_filename_1> = Input filename with the APD's from the first cell")
        print("<input_filename_2> = Input filename with the APD's from the second cell")
        print("<BCL> = Basic cycle length (period)")
        print("-------------------------------------------------------------------------------------------------------------")
        return 1

    # Get user inputs
    apd_file_name_1 = sys.argv[1]
    apd_file_name_2 = sys.argv[2]
    bcl = sys.argv[3]

    # Get the transmembrane potential from the input file
    data_1 = np.genfromtxt(apd_file_name_1)
    data_2 = np.genfromtxt(apd_file_name_2)

    #print(data_1.ndim)
    #print(data_2.ndim)

    if (data_1.ndim != 2) or (data_2.ndim != 2):
        block = True
    else:
        block = False

    if (block == False):
        apds_cell_100 = data_1[:,0]
        maxv_cell_100 = data_1[:,1]

        apds_cell_400 = data_2[:,0]
        maxv_cell_400 = data_2[:,1]
    
        diff_apd_cell_100 = apds_cell_100[0] - apds_cell_100[1]
        diff_apd_cell_400 = apds_cell_400[0] - apds_cell_400[1]

        if (maxv_cell_400[0] < 0.0):
            print("[BCL=%s] Block propagation" % (bcl))
        elif (abs(diff_apd_cell_100) < 4.0):
            print("[BCL=%s] No alternans" % (bcl))
        elif ( (diff_apd_cell_100 > 0.0) and (diff_apd_cell_400 > 0.0) ) or ( (diff_apd_cell_100 < 0.0) and (diff_apd_cell_400 < 0.0) ):
            print("[BCL=%s] Concordant alternans" % (bcl))
        elif ( (diff_apd_cell_100 > 0.0) and (diff_apd_cell_400 < 0.0) ) or ( (diff_apd_cell_100 < 0.0) and (diff_apd_cell_400 > 0.0) ):
            print("[BCL=%s] Discordant alternans" % (bcl))
    else:
        print("[BCL=%s] Block propagation" % (bcl))


if __name__ == "__main__":
    main()
