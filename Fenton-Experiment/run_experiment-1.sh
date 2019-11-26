#!/bin/bash

run_prepacing () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Running PREPACING with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Prepacing 
    ./bin/FentonExperiment examples/experiment-1/prepacing:cable-$2_bcl-$1.ini 
    mv cable-$2_bcl-$1.sst steady_state/experiment-1/cable-$2_bcl-$1.sst 
    echo "****************************************************************************************************"
}

run_simulation () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Running SIMULATION with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Experiment
    ./bin/FentonExperiment examples/experiment-1/cable-$2_bcl-$1.ini
    echo "****************************************************************************************************"
}

run_apd_calculation () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Calculating APD with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Script
    python ./scripts/calc_apd_elnaz.py ./results/$2cm/$1ms/sv-$3.dat 2 $1 1 100 65 > ./results/$2cm/$1ms/apd-cell-$3.dat
    python ./scripts/calc_apd_elnaz.py ./results/$2cm/$1ms/sv-$4.dat 2 $1 1 100 65 > ./results/$2cm/$1ms/apd-cell-$4.dat
    echo "****************************************************************************************************"
}

run_propagation_velocity_calculation () {
    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $2cm] Calculating CV with $1ms BCL ..."
    echo "****************************************************************************************************"
    # Script
    python ./scripts/calc_propagation_velocity.py ./results $2 0.0025 $1
    echo "****************************************************************************************************"
}

run_alternans_type () {
    # Script
    python ./scripts/calc_alternans_type.py ./results/$2cm/$1ms/apd-cell-$3.dat ./results/$2cm/$1ms/apd-cell-$4.dat $1 >> ./results/$2cm/alternans-$2cm.dat
}

run_cable_simulation () {
    
    PACINGS=($(seq 280 -5 150))

    ./clear_results.sh

    for PACING in "${PACINGS[@]}"; do
        run_prepacing $PACING $1
    done

    for PACING in "${PACINGS[@]}"; do
        run_simulation $PACING $1
    done

    for PACING in "${PACINGS[@]}"; do
        run_apd_calculation $PACING $1 $2 $3
    done

    echo "****************************************************************************************************"
    echo "[CABLE_LENGTH = $CABLE_LENGTH cm] Calculating alternans type ..."
    echo "****************************************************************************************************"
    for PACING in "${PACINGS[@]}"; do
        run_alternans_type $PACING $1 $2 $3
    done
    echo "****************************************************************************************************"

    #mkdir ./scripts/$1cm
    #mv ./scripts/*ms ./scripts/$1cm/
    #mv ./scripts/alternans-$1cm.dat ./scripts/$1cm/

    for PACING in "${PACINGS[@]}"; do
        run_propagation_velocity_calculation $PACING $1 $2 $3
    done

}

CABLE_LENGTHS=("5" "4.5" "4" "3.5" "3" "2.5" "2" "1.5" "1")
OFFSET_100=("400" "360" "320" "280" "240" "200" "160" "120" "80")
OFFSET_400=("1600" "1440" "1280" "1120" "960" "800" "640" "480" "320")
NUMBER_TESTS="8"

#CABLE_LENGTHS=("5")
#OFFSET_100=("400")
#OFFSET_400=("1600")
#NUMBER_TESTS="0"

for i in $(seq 0 1 $NUMBER_TESTS); do
    CABLE_LENGTH=${CABLE_LENGTHS[$i]}
    OFFSET_1=${OFFSET_100[$i]}
    OFFSET_2=${OFFSET_400[$i]}

    run_cable_simulation $CABLE_LENGTH $OFFSET_1 $OFFSET_2
done

