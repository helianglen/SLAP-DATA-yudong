% this section pulls the M2 power data and timetags
% from the RADTELEM files and creates mat files with the data
clc
clear all
close all

filesRad = dir('RAD*.slapbin');


for i = 1:length(filesRad)

  if filesRad(i).bytes  > 0   % skip empty files 
    disp(['Processing ', filesRad(i).name])  
    tic
    SLAP_read_raw_radiometer_data(filesRad(i).name)
    toc
  end
end

quit

