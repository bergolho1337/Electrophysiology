# This script will copy the files of the 'steady_state' and 'example' folders. 
# In order to allow the user to reproduce the results from the  simulations.
# ----------------------------------------------------------------------------------
#!/bin/bash

# Remove the current '.sst' files from the 'steady_state' folder
echo "[!] Removing old '.sst' files from 'steady_state' folder ..."
rm ../../../steady_state/*.sst

# Copy the '.sst' files from the current test
echo "[!] Copying the '.sst' files from the 5cm 'steady_state' folder ..."
cp ./steady_state/*.sst ../../../steady_state

# Remove the current '.ini' files from the 'examples' folder
echo "[!] Removing old '.ini' files from 'examples' folder ..."
rm ../../../examples/*.ini

# Copy the '.ini' files from the current test
echo "[!] Copying the '.ini' files from the 5cm 'examples' folder ..."
cp ./examples/*.ini ../../../examples
