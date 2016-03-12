% this section pulls the M2 power data and timetags
% from the RADTELEM files and creates mat files with the data
clc
clear all
%close all

% data_path = '~/Downloads/20150514 Sky Cal Data'
%data_path = '~/Desktop/20151014T190316 Dulles - Buoy - RTB'
%data_path = '/Volumes/SLAP DATA/20151106T174214 SLAPex CAN Afternoon Flight'
data_path = '.' 

script_path = pwd;

%cd(data_path);
filesRad = dir('RAD*.slapbin');
%cd(script_path);

for i = 1:length(filesRad)
    SLAP_ingest_radiometer(fullfile(data_path, filesRad(i).name), 'fb')
end

%cd(data_path);
%filesRad = dir('RDRRETURN*.slapbin');
%cd(script_path);

%for i = 1:length(filesRad)
%    SLAP_ingest_radar(fullfile(data_path, filesRad(i).name))
%end

beep
pause(1)
beep
pause(1)
beep
pause(1)
beep
pause(1)
beep
