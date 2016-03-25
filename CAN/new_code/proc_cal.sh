#!/bin/bash

#  process calibration data only

# Hybrid code from Albert's and Hemanshu's
# The data processing flow. The script should be run under each flight data directory  
# Before running this, make sure the OXTS*.mat data corresonding to the current experiment
# is located in the same directory. 
# Example
# cd 20151106T174214
# ../new_code/proc.sh 

 n_fb=`ls RADTELEM_*_FB.mat |wc -l`
 if [ $n_fb -eq 0 ]; then 
  # RADTELEM_20151108T210000.slapbin -> RADTELEM_20151108T210000_FB.mat
  matlab -nodesktop -nosplash -r new_run_L1B_allfiles
 else 
  echo FB.mat data processed ... skip to next step 
 fi
 matlab -nodesktop -nosplash -r m1_m2_corrector


exit 

