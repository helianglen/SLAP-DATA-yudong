#!/bin/bash

# Hybrid code from Albert's and Hemanshu's
# The data processing flow. The script should be run under each flight data directory  
# Before running this, make sure the OXTS*.mat data corresonding to the current experiment
# is located in the same directory. 
# Example
# cd 20151106T174214
# ../new_code/proc.sh 

 n_fb=`ls RDRRETUR*.mat |wc -l`
 if [ $n_fb -eq 0 ]; then 
  matlab -nodesktop -nosplash -r new_run_L1B_allfiles_radar 
 else 
  echo radar.mat data processed ... skip to next step 
 fi

 matlab -nodesktop -nosplash -r combine_step

 matlab -nodesktop -nosplash -r radar_6deg_universal

exit 


