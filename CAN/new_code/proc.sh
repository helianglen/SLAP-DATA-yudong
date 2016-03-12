#!/bin/bash

# Hybrid code from Albert's and Hemanshu's
# The data processing flow. The script should be run under each flight data directory  
# Before running this, make sure the OXTS*.mat data corresonding to the current experiment
# is located in the same directory. 
# Example
# cd 20151106T174214
# ../new_code/proc.sh 

 matlab -nodesktop -nosplash -r new_run_L1B_allfiles
 matlab -nodesktop -nosplash -r m1_m2_corrector
 matlab -nodesktop -nosplash -r combineAz
 matlab -nodesktop -nosplash -r combineADAQ
 matlab -nodesktop -nosplash -r consolidate_radiometer
 matlab -nodesktop -nosplash -r plot_radiometer


exit 

#----------------------------------------------------------------





# below is th original, based on Hemanshu's code 


# matlab -nodesktop -nosplash -r new_run_L1B_allfiles
 matlab -nodesktop -nosplash -r m1_m2_corrector
# matlab -nodesktop -nosplash -r combineAz
# matlab -nodesktop -nosplash -r combineADAQ
 matlab -nodesktop -nosplash -r radiometer_data_consolidation_CAN
 matlab -nodesktop -nosplash -r radiometer_data_plotter_CAN

