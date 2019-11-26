# =================================================================================================================================================================
# Author: Lucas Berg 
#
# This script calculates the range of combinations (sigma_c,G_gap) which will give us a prescribed 'sigma_bar' by following equation (1):
#
# sigma_bar = ( l ) / ( (l / sigma_c) + (A / G) ) .................................. (1)
#
# =================================================================================================================================================================
import numpy as np
import matplotlib.pyplot as plt

# Global variables
l = 0.01                    # {cm}
A = 235.6 * 1.0e-08         # {cm^2} 

# Output will be given in {mS/cm}
def calc_sigma_bar (s,G):
    return ( l ) / ( (l / s) + (A / G) )

# Output will be given in {mS}
def calc_G_gap (s_bar,s):
    return (-s_bar*s*A)/(s_bar*l-l*s)

def main ():

    # Flavio's Fenton experiment conductivity value
    sigma_bar = 0.035           # {mS/cm}

    # Setting the citoplasmatic conductivity range
    min_sigma_c = 1             # {mS/cm}
    max_sigma_c = 10            # {mS/cm}

    sigma_c = np.linspace(min_sigma_c,max_sigma_c,100)

    G_gap = calc_G_gap(sigma_bar,sigma_c)

    print("%30s\t%30s\t%30s" % ("sigma_c{mS/cm}","G_gap{uS}","sigma_bar{mS/cm}"))
    for i in range(len(sigma_c)):
        print("%30.20lf\t%30.20lf\t%30.20lf" % (sigma_c[i],G_gap[i]*1.0e+03,calc_sigma_bar(sigma_c[i],G_gap[i])))
    
    plt.plot(sigma_c,G_gap)
    plt.show()


if __name__ == "__main__":
    main()
