#!/bin/bash

CABLE_LENGTHS=("5" "4.5" "4" "3.5" "3" "2.5" "2" "1.5" "1")
SIZE="8"
PACINGS=($(seq 280 -5 150))

for i in $(seq 0 1 $SIZE); do
    CABLE_LENGTH=${CABLE_LENGTHS[$i]}
    mkdir results/${CABLE_LENGTH}cm
    for PACING in "${PACINGS[@]}"; do
        mkdir results/${CABLE_LENGTH}cm/${PACING}ms
    done
done

