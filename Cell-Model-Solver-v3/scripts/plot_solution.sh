#!/bin/bash

if test "$#" -ne 1; then
    echo "============================================================="
    echo "[!] ERROR! Illegal number of parameters"
    echo "============================================================="
    echo "Usage:> ./plot_solution.sh <solution_directory>"
    echo "-------------------------------------------------------------"
    echo "Example: ./plot_solution.sh tmp_mitchell_shaeffer"
    echo "============================================================="
    exit 1
fi

FILE_PATH="$1/V_t" 

# Call the plot script
python plot_ap.py $FILE_PATH