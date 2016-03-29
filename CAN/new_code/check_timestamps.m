%% plot timestamps from all the clock sources 
clear

%sometimes there is block of timetags that are bad

timeRad_all = deal([]); 

filesRad = dir('RAD*FB_m2data.mat');

% into a matlab .mat file first!!!
oxfile = dir('OXTSP*.mat');
load(oxfile.name)
timeGeo = gpstime;

% 16 second offset (converted to fraction of day)
% between the OxTS which is in GPS time and the rest of the
% data which is in UTC time

%YDT: 17 after july 2015
timeGeo = timeGeo + 17/86400;

% Currently we're using a linear two point cal between the cold sky obs on April 21, 2014 and
% the foam box on the day of the flight.


% NOTE: get azimuth scan angle data from the encoder files because there is an offset between the
% encoder and resolver (which has data in MOTIONTELEM files)


load AzData % This workspace variable is created in combineAz.m 

%% iterate on number of radiometer observation ".mat" files to process

for i = 1:length(filesRad)
    % displays which file is currently being processed
    filesRad(i).name(1:end-11)
    % the M1 corrected M2 counts (h and v pol) need to be generated first.
    % Then the 1.09 Hz noise is removed with "noise_remover.m". Finally, each h
    % and v pol data set is loaded here from "...m2data.m" as h2new_v2 and v2new_v2.

    load(filesRad(i).name)

    timeRad_all = [timeRad_all; time_file]; 

end 
    
   plot(timeRad_all, 'o') 
   datetick('y', 15) 
   hold on
   plot(timeGeo, ':')
   plot(timeMotAll, 'Color', [1 0 0]) 


