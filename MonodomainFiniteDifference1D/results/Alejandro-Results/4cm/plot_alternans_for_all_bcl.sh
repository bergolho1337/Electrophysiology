#!/bin/bash

ms="ms"

for bcl in $(seq 280 -5 100); do
	echo "[!] Plotting BCL = $bcl ..."
	folder_path="$bcl$ms"	
	python plot_potential.py $folder_path/sv-0.dat $folder_path/sv-0.pdf
	python plot_potential.py $folder_path/sv-240.dat $folder_path/sv-320.pdf
done
