%%
% determine the flight day of month
flight_directory = pwd;
flight_day = str2double((flight_directory(10:11)));
flight_month = str2double((flight_directory(8:9)));

% set time period of interest and load GPS receiver data
if flight_day == 12
    timeRadInit = datenum('2/12/2015 17:28:34');
    timeRadFinal = datenum('2/12/2015 17:55:55');
elseif flight_day == 25
    timeRadInit = datenum('2/25/2015 21:40');
    timeRadFinal = datenum('2/25/2015 23:07');
elseif flight_day == 1 && flight_month == 4
%     timeRadInit = datenum('4/1/2015 20:00:51');
%     timeRadFinal = datenum('4/1/2015 20:17:02');
    timeRadInit = datenum('4/1/2015 19:08:47');
    timeRadFinal = datenum('4/1/2015 19:11');
end



%% load and filter radiometer data

% filesRAD = dir('RAD*m2data.mat');
% [time_total, h_total, v_total]  = deal([]);
% for i = 14:15%2:length(filesRAD)
%     
%     load(filesRAD(i).name)
%     
%     h_total = [h_total; h2ant];
%     v_total = [v_total; v2ant];
%     time_total = [time_total; time_file];
% end
% save radiometer_data *total
% % %%
% 
load radiometer_data
% separate data for time period of interest
indRad = find(time_total >= timeRadInit & time_total < timeRadFinal );
timeRad = time_total(indRad);
hRad = h_total(indRad);
vRad = v_total(indRad);

%% get GPS information for time period of interest
% pre-process GPS data in excel file using importOxTSPostprocessed_v2.m
% into a matlab .mat file which is loaded next
geo = dir('Oxts*.mat');
load(geo(1).name)
% do some pre-processing of the GPS receiver time tags
timeGeo = gpstime;

% 16 timeond offset (converted to fraction of day)
% between the OxTS which is in GPS time and the rest of the
% data which is in UTC time
timeGeo = timeGeo + 16/86400;

% only get data from time period of interest
timeGeoInd = find(timeGeo >= timeRadInit & timeGeo < timeRadFinal );
timeGeoSlice = timeGeo(timeGeoInd);
lat = PosLatdeg(timeGeoInd);
lon = PosLondeg(timeGeoInd);
alt = PosAltm(timeGeoInd);
roll = AngleRolldeg(timeGeoInd);
pitch = AnglePitchdeg(timeGeoInd);
hdg = AngleHeadingdeg(timeGeoInd);
trk = AngleTrackdeg(timeGeoInd);

%% get scan angle data for time period of interest
% the azimuth scan angles were measured pos CCW from due North so we
% need to subtract them from 360 degrees to get them measured pos CW like
% heading and track angles which follow aircraft norms.

% The scan angle data is imported in CombineAncillaryData.m.
load Resolver_AzData 
az = 360 - azAll;

indMot = find(timeMotAll >= timeRadInit & timeMotAll < timeRadFinal );
timeMot = timeMotAll(indMot);
az = az(indMot);
%% interpolate GPS data and azimuth scan angle data to Rad time tag frequency

% This for-loop changes az values to decrease constantly as
% opposed to recycling to zero at 360 degrees. This is for interpolation
% purposes. This is computationally expensive so it is advantageous to only
% do it once then save the data for future calculations. Make sure to
% remove the az_cts file if you are processing a different time period of data
% so a new az_cts file can be created.

% cycind = find(diff(az)>0)+1;
% tic
% offset = 0;
% for m = 1:length(az)
%     if isempty(find(m == cycind, 1))
%         az(m) = az(m) - offset;
%     elseif ~isempty(find(m == cycind, 1))
%         offset = offset + 360;
%         az(m) = az(m) - offset;
%     end
% end
% toc

% interpolate az angles to Radiometer time frequency using timeRad as new index
az_interp = interp1(timeMot,az,timeRad);

% now interpolate OxTS components using timeGeoSlice as the original
% index and timeRad as the interpolating index
alt_interp  = interp1(timeGeoSlice, alt, timeRad);
pitch_interp = interp1(timeGeoSlice, pitch, timeRad);
roll_interp  = interp1(timeGeoSlice, roll, timeRad);
trk_interp  = interp1(timeGeoSlice, trk, timeRad);
hdg_interp  = interp1(timeGeoSlice, hdg, timeRad);
lonexpanded = interp1(timeGeoSlice, lon, timeRad);
latexpanded = interp1(timeGeoSlice, lat, timeRad);

% get altitude above ground level based on the national elevation database
% (NED) elevation for the each lat/lon coordinate
ned_vals = findElevData(latexpanded,lonexpanded);
AGL_interp = alt_interp - ned_vals;

% use SLAP's elevation angle, which changes with roll angle, to determine surface radius, srad, from
% nadir to location of observation
elev = 40 + roll_interp;
surface_radius = tand(elev).*AGL_interp;

% calculate complete az angle from true north. The azimuth offset = 0 in the SLAP measurements.
az_total = az_interp + hdg_interp ;

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

%% find indices of az values that are within 6 degree increments
increment = 134;% ~= 6 deg in scan angle

azcnt = 0;
azind = cell(2,1);

numlines = ceil(length(az_total)/increment);
[azavg,latavg, lonavg, havg, vavg, altavg, trkavg, rollavg, hdgavg, timeavg] = deal(zeros(numlines-1,1));
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
        havg(h) = nanmean(hRad(index));
        vavg(h) = nanmean(vRad(index));
end

% scale angle data back into range of [0 360] deg
for m = 1:length(azavg)
    azavg(m) = mod(azavg(m),360);
end

% save data to be processed by radiometer_plotter.m script
save('radiometer_6deg_avg.mat', 'hdgavg', 'azavg', 'trkavg',  'rollavg', 'lonavg', 'latavg', 'altavg', 'havg', 'vavg', 'azind', 'timeavg')
