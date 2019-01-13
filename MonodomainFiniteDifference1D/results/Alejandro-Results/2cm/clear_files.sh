#!/bin/bash

ms="ms"

for bcl in $(seq 280 -5 100); do
	echo "[!] Removing files from BCL = $bcl ..."
	folder_path="$bcl$ms"
	rm $folder_path/sv-60.dat
	rm $folder_path/sv-100.dat
	rm $folder_path/sv-180.dat
	rm $folder_path/sv-200.dat
	rm $folder_path/sv-240.dat
	rm $folder_path/sv-300.dat
	rm $folder_path/sv-400.dat
	rm $folder_path/sv-320.dat
done
