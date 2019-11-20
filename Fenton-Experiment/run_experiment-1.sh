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
    python ./scripts/calc_apd_elnaz.py ./scripts/$1ms/sv-$3.dat 2 $1 1 100 65 > ./scripts/$1ms/apd-cell-$3.dat
    python ./scripts/calc_apd_elnaz.py ./scripts/$1ms/sv-$4.dat 2 $1 1 100 65 > ./scripts/$1ms/apd-cell-$4.dat
    echo "****************************************************************************************************"
}

run_propagation_velocity_calculation () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Calculating CV with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Script
    python ./scripts/calc_propagation_velocity.py ./scripts/lucas-results $2 0.01 $1
    echo "****************************************************************************************************"
}

run_alternans_type () {
    # Script
    python ./scripts/calc_alternans_type.py ./scripts/$1ms/apd-cell-$3.dat ./scripts/$1ms/apd-cell-$4.dat $1 >> ./scripts/alternans-$2cm.dat
}

run_cable_simulation () {
    
    PACINGS=($(seq 280 -5 150))

    #./clear_results.sh

    #for PACING in "${PACINGS[@]}"; do
    #    run_prepacing $PACING $1
    #done

    #for PACING in "${PACINGS[@]}"; do
    #    run_simulation $PACING $1
    #done

    #for PACING in "${PACINGS[@]}"; do
    #    run_apd_calculation $PACING $1 $2 $3
    #done

    #echo "****************************************************************************************************"
    #echo "[CABLE_LENGTH = $CABLE_LENGTH cm] Calculating alternans type ..."
    #echo "****************************************************************************************************"
    #for PACING in "${PACINGS[@]}"; do
    #    run_alternans_type $PACING $1 $2 $3
    #done
    #echo "****************************************************************************************************"

    for PACING in "${PACINGS[@]}"; do
        run_propagation_velocity_calculation $PACING $1 $2 $3
    done

    #mkdir ./scripts/$1cm
    #mv ./scripts/*ms ./scripts/$1cm/
    #mv ./scripts/alternans-$1cm.dat ./scripts/$1cm/

}

CABLE_LENGTHS=("5" "4.5" "4" "3.5" "3" "2.5" "2" "1.5" "1")
OFFSET_100=("100" "90" "80" "70" "60" "50" "40" "30" "20")
OFFSET_400=("400" "360" "320" "280" "240" "200" "160" "120" "80")
NUMBER_TESTS="8"

for i in $(seq 0 1 $NUMBER_TESTS); do
    CABLE_LENGTH=${CABLE_LENGTHS[$i]}
    OFFSET_1=${OFFSET_100[$i]}
    OFFSET_2=${OFFSET_400[$i]}

    run_cable_simulation $CABLE_LENGTH $OFFSET_1 $OFFSET_2
done

