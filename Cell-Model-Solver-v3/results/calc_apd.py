import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Build a tuple where: T{i} = (t_i,dvdt_i)
def calc_derivatives (t,vm):
    n = len(vm)-1
    dvdt = {}
    for i in range(n):
        dvdt[t[i]] = vm[i+1]-vm[i] 
    return dvdt

def calc_apd (sv):
    # Copy the state vector data to local variables
    t = sv[:,0]
    vm = sv[:,1]
    
    # Calculate APD resting value
    min_value = np.amin(vm)
    max_value = np.amax(vm)
    APD_90_value = abs(max_value - min_value)*0.1 + min_value

    # Calculate the derivatives for each interval and sort then in decreasing order
    dvdt = calc_derivatives(t,vm)
    max_dvdt = sorted(dvdt.iteritems(), key=lambda (k,v): (v,k), reverse=True)

    for i in range(10):
        print max_dvdt[i]

    # Get the peak times of each action potential
    #default: peak_times = max_dvdt[0:2]
    peak_times = max_dvdt[3:5]
    peak_times.sort()
    
    for i in range(len(peak_times)):
        start_index = 0
        end_index = 0
        t1 = peak_times[i][0]
        t2 = peak_times[i][0]

        # Find the time where this stimulus occur in the data
        for j in range(len(t)):
            diff = abs(t[j] - t1)
            if (diff < 1.0e-08):
                start_index = j
        # Until the transmembrane potential reaches the APD_90_value we iterate over the data
        # When the difference between the current potential and the APD_value is minimum and when
        # we do not pass the time of the next stimulus we can save the end time of the APD
        min_diff = max_value
        min_t = t1
        for j in range(start_index,len(t)):
            # If the current time pass the start time of the next stimulus we break the loop
            if (i < len(peak_times)-1 and t[j] > peak_times[i+1][0]):
                min_t = t[j]
                end_index = j
                break
            # Until then we find the minimum difference between the APD value and the potential
            diff = abs(APD_90_value - vm[j])
            # Save the time when the difference is minimum
            if (diff < min_diff):
                min_diff = diff
                end_index = j
                min_t = t[j]
        t2 = min_t
        
        print("APD %d" % i)
        print("t1 = %.10lf" % t1)
        print("t2 = %.10lf" % t2)
        print("start_index = %d" % start_index)
        print("end_index = %d" % end_index)
        print("value = %.10lf" % (t2-t1))
        print("-------------------------------------------------------")
            



def main():

    if (len(sys.argv) != 2):
        print("Usage:> %s <input_file>" % (sys.argv[0]))
        sys.exit(1)
    else:
        input_file = sys.argv[1]

        data = np.genfromtxt(input_file)
        
        calc_apd(data)        
    

if __name__ == "__main__":
    main()