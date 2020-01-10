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

    os.chdir("outputs")
    output_file = open("../output_table.dat","w")
    for filename in glob.glob("atpi_*"):
        file = open(filename,"r")

        atpi = float(file.readline())
        Ko = float(file.readline())
        Ki = float(file.readline())
        Vm_mod = float(file.readline())
        GNa_mod = float(file.readline())
        GCaL_mod = float(file.readline())
        INaCa_mod = float(file.readline())

        Vm = float(file.readline())
        M = float(file.readline())
        H = float(file.readline())
        J = float(file.readline())
        Xr1 = float(file.readline())
        Xs = float(file.readline())
        S = float(file.readline())
        F = float(file.readline())
        F2 = float(file.readline())
        D_inf = float(file.readline())
        R_inf = float(file.readline())
        Xr2_inf = float(file.readline())

	hypoxia_mod = float(file.readline())
	hyperkelemia_mod = float(file.readline())
	acidosis_mod = float(file.readline())

        file.close()

	output_file.write("%g %g %g %g %g %g %g %g %g %g\n" % (hypoxia_mod,hyperkelemia_mod,acidosis_mod,atpi,Ko,Ki,Vm_mod,GNa_mod,GCaL_mod,INaCa_mod))

    output_file.close()

if __name__ == "__main__":
        main()
