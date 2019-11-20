#!/bin/bash

run_steady_state_simulation () {
    echo "****************************************************************************************************"
    echo "[!] Running steady-state simulation with (cable_length=$1 and sigma=$2) ..."
    echo "****************************************************************************************************"    
    ./bin/FentonExperiment examples/sst-cable-$1um-sigma-$2.ini
    cp output.sst steady_state/sst-cable-$1um-sigma-$2.sst

}

run_simulation () {
    echo "****************************************************************************************************"
    echo "[!] Running simulation with (cable_length=$1 and sigma=$2) ..."
    echo "****************************************************************************************************"    
    ./bin/FentonExperiment examples/cable-$1um-sigma-$2.ini
    mkdir results/experiment-1/cable-$1um-sigma-$2
    cp output/*.dat results/experiment-1/cable-$1um-sigma-$2
    ./clear_outputs.sh

}

SIGMAS=( 0.000100 0.000200 0.000300 0.000400 0.000500 0.000600 0.000700 0.000800 0.000900 )
CABLE_LENGTHS=( 15000.0 20000.0 25000.0 30000.0 35000.0 40000.0 45000.0 50000.0 )

#for SIGMA in "${SIGMAS[@]}"; do
#    for CABLE_LENGTH in "${CABLE_LENGTHS[@]}"; do
#	run_steady_state_simulation $CABLE_LENGTH $SIGMA
#    done
#done

for SIGMA in "${SIGMAS[@]}"; do
    for CABLE_LENGTH in "${CABLE_LENGTHS[@]}"; do
	run_simulation $CABLE_LENGTH $SIGMA
    done
done
