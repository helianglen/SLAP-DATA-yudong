% this section pulls the M2 power data and timetags
% from the RADTELEM files and creates mat files with the data
clc
clear all
%close all

data_path = '.' 

filesRad = dir('RDRRETURN*.slapbin');

for i = 1:length(filesRad)
    SLAP_ingest_radar(fullfile(data_path, filesRad(i).name))
end

quit

