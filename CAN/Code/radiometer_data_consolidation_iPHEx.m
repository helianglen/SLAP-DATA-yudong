%%
clear
filesRad = dir('RAD*FB_m2data.mat');

% determine flight date
flight_directory = pwd;
flight_day = str2double((flight_directory(10:11)));
flight_month = str2double((flight_directory(8:9)));
flight_year = str2double((flight_directory(4:7)));

files_to_process = 1:length(filesRad);

sample_period = 500e-6; %sec
sample_period_days = sample_period/86400;

% specify what data set is current being worked on to append to filenames
% of output
name = '_v4';

%%

% must pre-process GPS data in excel file using importOxTSPostprocessed_v2.m
% into a matlab .mat file first!!!
oxfile = dir('OxTSP*');
load(oxfile.name)
timeGeo = time;

% 16 second offset (converted to fraction of day)
% between the OxTS which is in GPS time and the rest of the
% data which is in UTC time
timeGeo = timeGeo + 16/86400;

% Currently we're using a linear two point cal between the cold sky obs on April 21, 2014 and
% the foam box on the day of the flight.


% NOTE: get azimuth scan angle data from the encoder files because there is an offset between the
% encoder and resolver (which has data in MOTIONTELEM files)
load AzData % This workspace variable is created in CombineAncillaryData.m
    
% the azimuth scan angles were measured pos CCW from due North so we
% need to subtract them from 360 degrees to get them measured pos CW like
% heading and track angles which follow the aircraft campaign norms
az = 360 - azAll;

%% iterate on number of radiometer observation ".mat" files to process

for i = files_to_process
    % displays which file is currently being processed
    filesRad(i).name(1:end-11)
    % the M1 corrected M2 counts (h and v pol) need to be generated first.
    % Then the 1.09 Hz noise is removed with "noise_remover.m". Finally, each h
    % and v pol data set is loaded here from "...m2data.m" as h2new_v2 and v2new_v2.

    load(filesRad(i).name(1:end-11), 'time')
    clear h2new_v2 v2new_v2
    load(filesRad(i).name)
    if exist('h2new_v2')
        h2ant = h2new_v2;
        v2ant = v2new_v2;
    else
        disp('Please generate 1.09 Hz noise removed data')
        break
    end
    
    timeRad = time;
    
    % find initial Rad timetag in sec
    flag_nonnan_first = find(~isnan(timeRad),1);
    timeRadInit = timeRad(flag_nonnan_first,:);
    
    % find final Rad timetag in sec
    flag_nonnan_final = find(~isnan(timeRad),1, 'last');
    timeRadFinal = timeRad(flag_nonnan_final,:);
    
    % % only use time tags of non-NaN data
    timeRad = timeRad(flag_nonnan_first:flag_nonnan_final);
    h2ant = h2ant(flag_nonnan_first*4-3:flag_nonnan_final*4);
    v2ant = v2ant(flag_nonnan_first*4-3:flag_nonnan_final*4);
    
    % find GPS data for time period of current ten minute radiometer data set
    timeGeoInd = find(timeGeo >= timeRadInit & timeGeo < timeRadFinal );
    timeGeoSlice = timeGeo(timeGeoInd);
    
    % get geolocation data at specified time range of current file
    lat = PosLatdeg(timeGeoInd);
    lon = PosLondeg(timeGeoInd);
    alt = PosAltm(timeGeoInd);
    roll = AngleRolldeg(timeGeoInd);
    pitch = AnglePitchdeg(timeGeoInd);
    hdg = AngleHeadingdeg(timeGeoInd);
    trk = AngleTrackdeg(timeGeoInd);
   
    % The timetags in the radiometer file  have some peculiar traits. They
    % jump ephemerally every once in a while, both to lower
    % and higher values. Additionally, they jump up permanently every once in a
    % while so the final time is larger in value than what would be
    % expected if the time was steadily increasing every 2 ms as it should be.
    
    % call to function to "repair" time tags from 2014 flight data
    time = repair_time_tags(timeRad);
    
    % expand the timetag for the RAD file to account for the 4 PRIs at
    % timetag + 0, timetag + 0.5 ms, timetag + 1 ms, and timetag + 1.5 ms
    time_file = [time, time + sample_period_days, time + sample_period_days*2, time + sample_period_days*3];
    timeRad =  reshape(time_file', numel(time_file), 1);
    
    indRad = find(timeRad >= timeMotAll(1) & timeRad < timeMotAll(end) );
    timeRad = timeRad(indRad);
    h2ant = h2ant(indRad);
    v2ant = v2ant(indRad);
    
    % reset az values less than zero to be within [0 360]
    az(az<0) = az(az<0) + 360;
    
    % interpolate az angles to Rad time frequency using timeRad as new index
    az_interp = interp1(timeMotAll,az,timeRad);
    
    % get altitude about ground level based on the national elevation database
    % (NED) elevation for the each lat/lon coordinate
    ned_vals = findElevData(lat,lon);
    AGL = alt - ned_vals;
    
    % now interpolate OxTS components using timeGeoSlice as the original
    % index and timeRad as the interpolating index
    alt_interp = interp1(timeGeoSlice, AGL, timeRad);
    pitch_interp = interp1(timeGeoSlice, pitch, timeRad);
    roll_interp = interp1(timeGeoSlice, roll, timeRad);
    trk_interp = interp1(timeGeoSlice, trk, timeRad);
    hdg_interp = interp1(timeGeoSlice, hdg, timeRad);
    lonexpanded = interp1(timeGeoSlice, lon, timeRad);
    latexpanded = interp1(timeGeoSlice, lat, timeRad);

    % use SLAP's incidence angle, which changes with roll angle, to determine surface radius, srad, from
    % nadir to location of observation
    elev = 40 - roll_interp .* sind(az_interp);
    surface_radius=tand(elev).*alt_interp;
    range=alt_interp./cosd(elev);
    
    % calculate complete az angle from true north. The azimuth offset = 0 in the SLAP measurements.
    az_total = az_interp + hdg_interp;
    
    % radius of earth in meters
    r_earth = 6371e3;
    
    % determine x and y offset from nadir (directly below the aircraft)
    % of each observation based on azimuth angle of SLAP at time of
    % observation.
    
    yoff = surface_radius .* cosd(az_total); % meters
    xoff = surface_radius .* sind(az_total); % meters
    
    % convert x and y offset into lat-lon offset, then add to current lat-lon coordinates.
    laty = yoff  .* 360./(2*pi*r_earth) + latexpanded ;	%pixel absolute x location (longitude)
    % horizontal offset values in degrees depend on latitude on Earth
    lonx = xoff  .* 360./(2*pi*r_earth.*cosd(latexpanded )) + lonexpanded ;	%pixel absolute y location (latitude)
    %% get 6 deg scan angle averages of radiometric, azimuthal scan angle and GPS location data
    
    increment = 134;% ~= 6 deg in scan angle
    
    azcnt = 0;
    azind = cell(2,1);
    
    numlines = ceil(length(az_total)/increment);
    [azavg,latavg, lonavg, havg, vavg, altavg, trkavg, rollavg, hdgavg, timeavg] = deal(zeros(numlines,1));
    for h = 1:numlines-1
        index = (increment*(h-1)+1):increment*h;
        azavg(h) = az_total(index(floor(end/2)));
        latavg(h) = laty(index(floor(end/2)));
        lonavg(h) = lonx(index(floor(end/2)));
        altavg(h) = alt_interp(index(floor(end/2)));
        trkavg(h) = trk_interp(index(floor(end/2)));
        hdgavg(h) = hdg_interp(index(floor(end/2)));
        rollavg(h) = roll_interp(index(floor(end/2)));
        timeavg(h) = timeRad(index(floor(end/2)));
        havg(h) = nanmean(h2ant(index));
        vavg(h) = nanmean(v2ant(index));
    end
    
    save([filesRad(i).name(1:end-11), name, '.mat'], 'hdgavg', 'azavg', 'trkavg',  'rollavg', 'lonavg', 'latavg', 'altavg', 'havg', 'vavg', 'timeavg')
end
