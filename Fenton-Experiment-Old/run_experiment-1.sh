#!/bin/bash

run_prepacing () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Running PREPACING with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Prepacing 
    ./bin/FentonExperiment examples/experiment-1/prepacing:cable-$2_bcl-$1.ini
    cp output.sst steady_state/experiment-1/cable-$2_bcl-$1.sst 
    #./clear_results.sh
    echo "****************************************************************************************************"
}

run_simulation () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Running SIMULATION with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Experiment
    ./bin/FentonExperiment examples/experiment-1/cable-$2_bcl-$1.ini
    mkdir scripts/$1ms
    cp output/*.dat scripts/$1ms
    echo "****************************************************************************************************"
}

run_apd_calculation () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Calculating APD with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Script
    python ./scripts/calc_apd_elnaz.py ./scripts/$1ms/sv-100.dat 2 $1 1 100 65 > ./scripts/$1ms/apd-cell-100.dat
    python ./scripts/calc_apd_elnaz.py ./scripts/$1ms/sv-400.dat 2 $1 1 100 65 > ./scripts/$1ms/apd-cell-400.dat
    echo "****************************************************************************************************"
}

run_alternans_type () {
    # Script
    python ./scripts/calc_alternans_type.py ./scripts/$1ms/apd-cell-100.dat ./scripts/$1ms/apd-cell-400.dat $1 >> ./scripts/alternans-$2cm.dat
}

#PACINGS=( 200 210 220 230 240 250 260 270 280 290 300 )
PACINGS=($(seq 280 -5 270))
CABLE_LENGTH="5"

./clear_results.sh

for PACING in "${PACINGS[@]}"; do
    run_prepacing $PACING $CABLE_LENGTH
done

for PACING in "${PACINGS[@]}"; do
    run_simulation $PACING $CABLE_LENGTH
done

for PACING in "${PACINGS[@]}"; do
    run_apd_calculation $PACING $CABLE_LENGTH
done

echo "****************************************************************************************************"
echo "[CABLE_LENGTH = $CABLE_LENGTH cm] Calculating alternans type ..."
echo "****************************************************************************************************"
for PACING in "${PACINGS[@]}"; do
    run_alternans_type $PACING $CABLE_LENGTH 
done
echo "****************************************************************************************************"