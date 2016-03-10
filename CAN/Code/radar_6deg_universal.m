%%
clear

% get flight date
flight_directory = pwd;
flight_day = str2double((flight_directory(10:11)));
flight_month = str2double((flight_directory(8:9)));
flight_year = str2double((flight_directory(4:7)));

%% Initialize key options
flagHalf = 1;
% flag to process data for fore or aft half-scan
fore = 1;

% flag to plot h or v pol data on Google Earth
GoogleEarth  = 1;
pol = 0; % 1 == hpol and 0 == vpol
co_pol = 0;

files = dir('RDRRET*00.mat');
files_to_process = 3:length(files);
%%
sample_period = 500e-6; %sec
sample_period_days = sample_period/86400;

if flight_year == 2014
    power_dBm = [-100 -95:2:-79];
    cal_hpol_counts = [300 800 1225 1900 3000 4600 7225 11450 18000 27375];
    cal_vpol_counts = [230 600 910 1410 2200 3435 5375 8390 13250 21300];
    timeRadInit = datenum('5/21/2014 18:23:44');
    timeRadFinal = datenum('5/21/2014 19:04:55');
    
elseif flight_year == 2015 && flight_month == 2
    power_dBm = [-85:5:-50, -49:-45, -40:-30];
    cal_hpol_counts = [1;3.2;10.7;33.5;106;335;1062;3342;4190;5270;6620;8270;10470;31700;38500;46700];
    % don't have an accurate radar calibration for V-pol for these flights
    % so just use H-pol calibration.
    cal_vpol_counts = cal_hpol_counts;
    if flight_day == 25
        timeRadInit = datenum('2/25/2015 22:25:58');
        timeRadFinal = datenum('2/25/2015 23:06:25');
    elseif flight_day == 12
        
    end
elseif flight_year == 2015 && flight_month >= 3
    % use SLAP Radar ground calibration data points from documentation
    power_dBm = [-85:5:-45, -43, -41, -40:-38]';
    
    cal_hpol_counts = [1; 3.1; 10.3; 32.5; 102.5; 322; 1017; 3190; 9980; 15600; 24100; 29800; 36500; 44000];
    cal_vpol_counts = [1;3.80000000000000;11.7000000000000;36.7000000000000;116.500000000000;366;1152;3640;11350;17750;27230;33450;40500;48300];
    timeRadInit = datenum('4/1/2015 20:31:18');
    timeRadFinal = datenum('4/1/2015 21:20:15');
end

power_mW = 10.^(power_dBm/10);

% only need power in dBm values before the radar return saturates. The
% cal_hpol_counts are only given before saturation so the number of valid data
% points of power_mW corresponds to the first n (= length(cal_hpol_counts) data points
power_mW_hpol = power_mW(1:length(cal_hpol_counts));

% same goes for v-pol
power_mW_vpol = power_mW(1:length(cal_vpol_counts));


noiseFloor_h = cal_hpol_counts(1);
noiseFloor_v = cal_vpol_counts(1);

saturationLimit_h = cal_hpol_counts(end);
saturationLimit_v = cal_vpol_counts(end);

% noise Floor and saturation Limit are intrinsically used since the
% interpolated method used to convert from raw counts to mW is a piecewise
% linear fit which turns raw count values outside the calibration range
% into NaN power values.

rangeGateDistance = 300; % meters

% initialize variables
f = 1254000000; % Hz
c = 299792458; % m/s
Pt = 10^(47/10);
Gt = 17.7^(47/10);
Gr = Gt;

alpha = 1.33;
phi = 40;

tau = 1e-6;
theta_ah = 13.4;
theta_av = 13.1;

% beam width information by polarization from SLAP documentation
beam_width_h_pol = 20.5;
beam_width_v_pol = 19.5;

% load azimuthal scan angle data
load AzData
azAll = 360 - azAll;

indMot = find(timeMotAll >= timeRadInit & timeMotAll < timeRadFinal );
timeMot = timeMotAll(indMot);
az = azAll(indMot);

% for the 2015 flights, load step attenuator data
if flight_year == 2015
    load step_all
end

%% only use GPS receiver data for time period of interest

% must pre-process GPS data in excel file using importOxTSPostprocessed_v2.m
% into a matlab .mat file first!!!
oxfile = dir('OxTS*.mat');
load(oxfile(1).name)

% 16 second offset (converted to fraction of day)
% between the OxTS which is in GPS time and the rest of the
% data which is in UTC time
timeGeo = time + 16/86400;

timeGeoInd = find(timeGeo >= timeRadInit & timeGeo < timeRadFinal );
timeGeoSlice = timeGeo(timeGeoInd);
lat = PosLatdeg(timeGeoInd);
lon = PosLondeg(timeGeoInd);
alt = PosAltm(timeGeoInd); % meters
roll = AngleRolldeg(timeGeoInd);
pitch = AnglePitchdeg(timeGeoInd);
hdg = AngleHeadingdeg(timeGeoInd);


%% precalculate altitude above ground level (AGL)
ned_vals = findElevData(lat,lon);
AGL = alt - ned_vals;

%% roll flag and limits
flagRoll = 0;
% Aircraft roll angle limits if removing data during a turn
rollmin = -5;
rollmax = 5;

%%
for k = files_to_process
    timestr = files(k).name(end-9:end-6)
    
    % load time, radar_h, and radar_v data
    load(files(k).name)
    timeRad = time;
    
    %% call to function to "repair" time tags
    
    if flight_year == 2014
        % call to function to "repair" time tags from 2014 radar flight data
        timeRad = repair_time_tags(timeRad);
        % create step attenuator data set of zeros for future calculations
        [step_h_pol, step_v_pol] = deal(zeros(length(timeRad),1));
    else
        % For 2015 data, the time tags are much better. However, they need to be
        % recreated still because of hardware issues that cause the time tags to
        % not continuously increase every once in a while. The instrument is
        % still scanning so this is a valid fix.
        time_new = time(1):0.0005/86400:(time(1)+((length(time)-1)*0.0005/86400));
        timeRad = time_new';
                
        % interpolate step attenuator data for 2015 flight data using previous neighbor method
        step_h_pol = interp1(step_all(:,1),step_all(:,2),timeRad, 'previous');
        step_v_pol = interp1(step_all(:,1),step_all(:,3),timeRad, 'previous');
    end
    %% interpolate GPS data and azimuth scan angle data to Rad time tag frequency
    
    % interpolate az angles to Rad time frequency using timeRad as new index
    az_interp = interp1(timeMot,az,timeRad);
    
    % now interpolate OxTS components using timeGeoSlice as the original
    % index and timeRad as the interpolating index
    AGL_interp = interp1(timeGeoSlice, AGL, timeRad);
    pitch_interp = interp1(timeGeoSlice, pitch, timeRad);
    roll_interp = interp1(timeGeoSlice, roll, timeRad);
    trk_interp = interp1(timeGeoSlice, trk, timeRad);
    hdg_interp = interp1(timeGeoSlice, hdg, timeRad);
    lonexpanded = interp1(timeGeoSlice, lon, timeRad);
    latexpanded = interp1(timeGeoSlice, lat, timeRad);
    
    % use SLAP's incidence angle, which changes with roll angle, to determine
    % the elevation angle
    elev_angle = phi - roll_interp .* sind(az_interp);
    
    % use SLAP's given elevation angle to determine surface radius, srad, from
    % nadir to location of observation
    surface_radius=tand(elev_angle).*AGL_interp;
    
    % precalculate slant range
    slant = AGL_interp./ cosd(elev_angle);
    
    % calculate complete az angle from true north. The azimuth offset = 0 in the SLAP measurements.
    az_total = az_interp + hdg_interp;% + az_offset;
    
    % radius of earth in meters
    r_earth = 6371e3;
    
    % determine x and y offset from nadir (directly below the aircraft)
    % of each observation based on azimuth angle of SLAP at time of
    % observation.
    
    yoff = surface_radius .* cosd(az_total); % meters
    xoff = surface_radius .* sind(az_total); % meters
    
    % convert x and y offset into lat-lon offset, then add to current lat-lon coordinates.
    lat_interp = yoff  .* 360./(2*pi*r_earth) + latexpanded ;	%pixel absolute x location (longitude)
    % horizontal offset values in degrees depend on latitude on Earth
    lon_interp = xoff  .* 360./(2*pi*r_earth.*cosd(latexpanded )) + lonexpanded ;	%pixel absolute y location (latitude)
    
    %% get sigma0 data for h and v pol at six degree scan angle averages
    
    % determine which data corresponds to co-pol and which corresponds to
    % cross-pol using the PCM vector.
    
    pcm_sep = bitget(pcm,3);
    % if the PCM has a 0 in its 4s bit, it's H transmitting.
    % if it has a 1 in its 4s bit, it's V transmitting.
    % the nomenclature is as follows: indstart_xy where x is the transmit pol
    % and y is the receive pol
    if pcm(1) == 0
        indstart_hh = 1;
        indstart_vh = 2;
        indstart_hv = 1;
        indstart_vv = 2;
    else
        indstart_hh = 2;
        indstart_vh = 1;
        indstart_hv = 2;
        indstart_vv = 1;
    end
    
    increment = 134;% ~= 6 deg in scan angle
    
    numlines = floor(length(az_total)/increment);
    
    % initialize variables
    [azavg, latavg, lonavg, mean_hh, mean_vv, mean_hv, mean_vh, altavg, trkavg, rollavg, hdgavg, timeavg] = deal(zeros(numlines,1));
    
    for j = 1:numlines
        % get radar data points for each range bin
        ind = (increment*(j-1)+1):increment*j;
        data_hh = radar_h( ind(indstart_hh:2:end),:);
        data_vv = radar_v( ind(indstart_vv:2:end),:);
        
        data_step_hh = step_h_pol( ind(indstart_hh:2:end));
        data_step_vv = step_v_pol( ind(indstart_vv:2:end));
        
        data_hv = radar_v( ind(indstart_hv:2:end),:);
        data_vh = radar_h( ind(indstart_vh:2:end),:);
        
        data_step_hv = step_v_pol( ind(indstart_hv:2:end));
        data_step_vh = step_h_pol( ind(indstart_vh:2:end));
        
        % use middle of current 6 degree data set indices for
        % geolocation and scan angle values
        indmid = ind(floor(length(ind)/2));
        
        latavg(j) = lat_interp(indmid);
        lonavg(j) = lon_interp(indmid);
        altavg(j) = AGL_interp(indmid);
        trkavg(j) = trk_interp(indmid);
        hdgavg(j) = hdg_interp(indmid);
        rollavg(j) = roll_interp(indmid);
        timeavg(j) = timeRad(indmid);
        azavg(j) = az_total(indmid);
        
        % initialize sigma0 values for 4 pol combinations for this 6 deg avg data set (=66 data points)
        [sigma0hh, sigma0vv, sigma0hv, sigma0vh] = deal(zeros(size(data_vv,1),1));
        
        % initialize  radar footprint for H-pol
        footprintm = (slant(indmid).*tand(beam_width_h_pol/2))./cosd(elev_angle(indmid));
        area_h = pi .* footprintm^2;
        
        % initialize  radar footprint for V-pol
        footprintm = (slant(indmid)*tand(beam_width_v_pol/2))./cosd(elev_angle(indmid));
        area_v = pi .* footprintm^2;
        
        for i = 1:size(data_hh,1)
            %% This section calculates HH sigma0
            if pol && co_pol
                % determine indices of counts above the noise floor
                indCountsWithinBounds = find(data_hh(i,:)> noiseFloor_h);
                % if there is no radar return data above the noise floor, the
                % code skips data point and leaves sigma0h's default value of 0
                if ~isempty(indCountsWithinBounds)
                    % get data points above the noise floor
                    usefulCounts = data_hh(i,indCountsWithinBounds);
                    
                    % distance from range gate 1 to each range gate with usable
                    % radar return data
                    surface_ranges = (indCountsWithinBounds-1) * rangeGateDistance;
                    % determine total surface distance
                    surfaceDistance = surface_radius(indmid) + surface_ranges;
                    % range from SLAP (at altitude) to each range gate is determined by the
                    % Pythagorean thereom.
                    R = sqrt ( surfaceDistance.^ 2 + altavg(j)^2);
                    
                    % interpolate counts values to determine power in mW
                    powerReceived_h_pol = interp1(cal_hpol_counts, power_mW, usefulCounts);
                    % convert power from mW to dBm to be able to add step atten
                    powerReceived_h_pol = 10*log10(powerReceived_h_pol);
                    
                    % factor in step attenuator data
                    powerReceived_h_pol = powerReceived_h_pol + data_step_hh(i);
                    
                    % convert power from dBm to mW
                    powerReceived_h_pol = 10.^(powerReceived_h_pol/10);
                    
                    % equation in linear (mW) space to solve for sigma
                    sumPandR = nansum(powerReceived_h_pol.*R.^4);
                    sigma_hh = sumPandR / ( Pt * Gt * ( c^2 / (4*pi*f)^2 ) * 1/(4*pi) * Gr);
                    
                    % sum non-nan values of (sigma/area) to get sigma0 in mW
                    sigma0hh(i) =  10*log10(  nansum(sigma_hh./area_h) );
                end
                %% This section calculates HV sigma0
            elseif pol && ~co_pol
                % determine indices of counts above the noise floor
                indCountsWithinBounds = find(data_hv(i,:)> noiseFloor_h);
                % if there is no radar return data above the noise floor, the
                % code skips data point and leaves sigma0h's default value of 0
                if ~isempty(indCountsWithinBounds)
                    % get data points above the noise floor
                    usefulCounts = data_hh(i,indCountsWithinBounds);
                    
                    % distance from range gate 1 to each range gate with usable
                    % radar return data
                    surface_ranges = (indCountsWithinBounds-1) * rangeGateDistance;
                    % determine total surface distance
                    surfaceDistance = surface_radius(indmid) + surface_ranges;
                    % range from SLAP (at altitude) to each range gate is determined by the
                    % Pythagorean thereom.
                    R = sqrt ( surfaceDistance.^ 2 + altavg(j)^2);
                    
                    % interpolate counts values to receive pol to determine power in mW
                    powerReceived_v_pol = interp1(cal_vpol_counts, power_mW, usefulCounts);
                    
                    % convert power from mW to dBm to be able to add step atten
                    powerReceived_v_pol = 10*log10(powerReceived_v_pol);
                    
                    % factor in step attenuator data
                    powerReceived_v_pol = powerReceived_v_pol + data_step_hv(i);
                    
                    % convert power from dBm to mW
                    powerReceived_v_pol = 10.^(powerReceived_v_pol/10);
                    % equation in linear (mW) space to solve for sigma
                    sumPandR = sum(powerReceived_v_pol.*R.^4);
                    sigma_hv = sumPandR / ( Pt * Gt * ( c^2 / (4*pi*f)^2 ) * 1/(4*pi) * Gr);
                    
                    % sum non-nan values of (sigma/area) to get sigma0 in mW.
                    % must divide by receiving pol footprint area.
                    sigma0hv(i) =  10*log10( nansum(sigma_hv./area_v));
                end
                %% This section calculates VV sigma0
            elseif ~pol && co_pol
                % determine indices of counts above the noise floor
                indCountsWithinBounds = find(data_vv(i,:)> noiseFloor_v);
                % if there is no radar return data above the noise floor, the
                % code skips data point and leaves sigma0vv's default value of 0
                if ~isempty(indCountsWithinBounds)
                    % get data points above the noise floor
                    usefulCounts = data_vv(i,indCountsWithinBounds);
                    
                    % distance from range gate 1 to each range gate with usable
                    % radar return data
                    surface_ranges = (indCountsWithinBounds-1) * rangeGateDistance;
                    % add these two to determine total surface distance
                    surfaceDistance = surface_radius(indmid) + surface_ranges;
                    % range from SLAP to each range gate is determined by the
                    % Pythagorean thereom.
                    R = sqrt ( surfaceDistance.^ 2 + altavg(j)^2);
                    
                    % interpolate counts values to determine power in mW\
                    powerReceived_v_pol = interp1(cal_vpol_counts, power_mW, usefulCounts);
                    
                    % convert power from mW to dBm to be able to add step atten
                    powerReceived_v_pol = 10*log10(powerReceived_v_pol);
                    % factor in step attenuator data
                    powerReceived_v_pol = powerReceived_v_pol + data_step_vv(i);
                    
                    % convert power from dBm to mW
                    powerReceived_v_pol = 10.^(powerReceived_v_pol/10);
                    % equation in linear (mW) space to solve for sigma
                    sumPandR = nansum(powerReceived_v_pol.*R.^4);
                    sigma_vv = sumPandR / ( Pt * Gt * ( c^2 / (4*pi*f)^2 ) * 1/(4*pi) * Gr);
                    
                    % sum non-nan values of (sigma/area) to get sigma0 in mW
                    sigma0vv(i) =  10*log10( nansum(sigma_vv./area_v));
                end
                %% This section calculates VH sigma0
            elseif ~pol && ~co_pol
                % determine indices of counts above the noise floor
                indCountsWithinBounds = find(data_vh(i,:)> noiseFloor_v);
                % if there is no radar return data above the noise floor, the
                % code skips data point and leaves sigma0vv's default value of 0
                if ~isempty(indCountsWithinBounds)
                    % get data points above the noise floor
                    usefulCounts = data_vh(i,indCountsWithinBounds);
                    
                    % distance from range gate 1 to each range gate with usable
                    % radar return data
                    surface_ranges = (indCountsWithinBounds-1) * rangeGateDistance;
                    % add these two to determine total surface distance
                    surfaceDistance = surface_radius(indmid) + surface_ranges;
                    % range from SLAP to each range gate is determined by the
                    % Pythagorean thereom.
                    R = sqrt ( surfaceDistance.^ 2 + altavg(j)^2);
                    
                    % interpolate counts values to receive pol to determine power in mW
                    
                    powerReceived_h_pol = interp1(cal_hpol_counts, power_mW, usefulCounts);
                    
                    % convert power from mW to dBm to be able to add step atten
                    powerReceived_h_pol = 10*log10(powerReceived_h_pol);
                    % factor in step attenuator data
                    powerReceived_h_pol = powerReceived_h_pol + data_step_vh(i);
                    
                    % convert power from dBm to mW
                    powerReceived_h_pol = 10.^(powerReceived_h_pol/10);
                    
                    % equation in linear (mW) space to solve for sigma
                    sumPandR = nansum(powerReceived_h_pol.*R.^4);
                    sigma_vh = sumPandR / ( Pt * Gt * ( c^2 / (4*pi*f)^2 ) * 1/(4*pi) * Gr);
                    
                    % sum non-nan values of (sigma/area) to get sigma0 in mW
                    % must divide by receiving pol footprint area.
                    sigma0vh(i) = 10*log10( nansum(sigma_vh./area_h));
                end
            end
        end
        % get six degree scan angle average of sigma0 for h and v pol
        mean_hh(j) = nanmean(sigma0hh);
        mean_hv(j) = nanmean(sigma0hv);
        mean_vh(j) = nanmean(sigma0vh);
        mean_vv(j) = nanmean(sigma0vv);
    end
    
    if co_pol
        copolstr = '_co_pol';
    else
        copolstr = '_cross_pol';
    end
   
    
    % save data for future processing
    %     save([files(k).name(1:end-4), '_6deg.mat'], 'mean_hh', 'mean_vv', 'latavg', 'lonavg', 'altavg', ...
    %         'mean_hv', 'mean_vh','hdgavg', 'rollavg', 'trkavg', 'timeavg', 'azavg')
    
    if GoogleEarth
       % only plot on either the fore or aft half scan if flagHalf is on      
        if flagHalf
            % initialize flag_half_scan as zeros
            flag_half_scan = zeros(length(trkavg),1);
            
            for k = 1:length(trkavg)
                % get 180 degrees of scan width around track angle to get front half
                scan_lower = trkavg(k) - 90;
                scan_upper = trkavg(k) + 90;
                scan_upper(scan_upper>360) = scan_upper(scan_upper>360) - 360;
                scan_lower(scan_lower<0) = scan_lower(scan_lower<0) + 360;
                
                % get indices of azimuth angles within the 180 degrees of the front
                % half
                flagf1 = azavg(k)<scan_upper;
                flagf2 = azavg(k)>=scan_lower;
                
                if trkavg(k) < 180;
                    if scan_lower > 180 && (flagf1 || flagf2)
                        flag_half_scan(k) = k;
                    elseif scan_lower <= 180 && (flagf1 && flagf2)
                        flag_half_scan(k) = k;
                    end
                else
                    if scan_upper > 180 && (flagf1 && flagf2)
                        flag_half_scan(k) = k;
                    elseif scan_upper <= 180 && (flagf1 || flagf2)
                        flag_half_scan(k) = k;
                    end
                end
            end
            
            % turn flag_half_scan into a vector of logicals where ones are indices of
            % fore half scan data points
            flag_half_scan = flag_half_scan>0;
            
            if ~fore
                % get indices of aft half scan as the opposite of the indices of the fore
                % half scan if the flag for the fore half scan is turned off
                flag_half_scan = ~flag_half_scan;
            end
            
        else
            % if flagHalf is off, to plot fore and aft half scans together, get all logical ones
            flag_half_scan = find(latavg>0);
        end
        
        % get geolocation data for fore or aft half scans using flag_half_scan
        trk_half_scan = trkavg(flag_half_scan);
        lon_half_scan = lonavg(flag_half_scan);
        lat_half_scan = latavg(flag_half_scan);
        roll_half_scan = rollavg(flag_half_scan);
        vv_half_scan = mean_vv(flag_half_scan);
        hh_half_scan = mean_hh(flag_half_scan);
        hv_half_scan = mean_hv(flag_half_scan);
        vh_half_scan = mean_vh(flag_half_scan);
        alt_half_scan = altavg(flag_half_scan);
        time_half_scan = timeavg(flag_half_scan);
        az_half_scan = azavg(flag_half_scan);
        hdg_half_scan = hdgavg(flag_half_scan);
        
        minBin= - 110;
        maxBin= -80;
        
        if pol
            beam_width = beam_width_h_pol;
            if co_pol
                poltype = '_HH';
                Tb = hh_half_scan;
            else
                poltype = '_HV';
                Tb = hv_half_scan;
            end
        else
            beam_width = beam_width_v_pol;
            if co_pol
                poltype = '_VV';
                Tb = vv_half_scan;
            else
                poltype = '_VH';
                Tb = vh_half_scan;
            end
        end
        
        Lat = lat_half_scan;
        Long = lon_half_scan;
        Alt = alt_half_scan;
        
        % if flagRoll is set on, perform the roll filter
        if flagRoll
            indRoll = find(roll_half_scan < rollmax & roll_half_scan > rollmin);
            Tb = Tb(indRoll);
            Lat = lat_half_scan(indRoll);
            Long = lon_half_scan(indRoll);
            Alt = alt_half_scan(indRoll);
        end
        
        % determine footprint size based on altitude
        slant = Alt./ cosd(phi);
        
        footprintm = (slant*tand(beam_width/2))./cosd(phi);
        
        % convert footprint radius in meters to degrees longitude
        footprint = footprintm / (pi*6378e3) * 180;
        
        % generate KML file to plot on Google Earth with lat/lon and Tb/SM data
        if fore
            name = '_Fore';
        else
            name = '_Aft';
        end
        
        fida=fopen(strcat(timestr, '_radar', poltype, name,  '.kml'),'w');
        
        line1=['<?xml version="1.0" encoding="UTF-8"?>'];
        line2=['<kml xmlns="http://earth.google.com/kml/2.0">'];
        line3=['<Document>'];
        line4=[strcat('<name>',strcat(timestr, '_radar', poltype, name),'</name>')];
        line6=['<Style id="air"><icon><href>root://icons/palette-4.png?x=160&amp;y=0&amp;w=32&amp;h=32</href></icon>'];
        line7=['<LineStyle><color>FF0000FF</color><width>2.0</width></LineStyle></Style>'];
        line8=['<Style id="gnd"><icon><href>root://icons/palette-4.png?x=160&amp;y=0&amp;w=32&amp;h=32</href></icon>'];
        line9=['<LineStyle><color>FFBBFF00</color><width>2.0</width></LineStyle></Style>'];
        line10=['<IconStyle><color>FFFFFFFF</color><scale>0.7</scale></IconStyle>'];
        line11=['<LabelStyle><color>FFFFFFFF</color><scale>0.7</scale></LabelStyle></Style>'];
        line_folder=['</Folder>'];
        line_last=['</Document></kml>'];
        
        cls=flipud(['7F0000AA';'7F0000D4';'7F0000FF';'7F002AFF';'7F0055FF';'7F0080FF';'7F00AAFF'; ...
            '7F00D4FF';'7F00FFFF';'7F2AFFD4';'7F54FFAA';'7F80FF80';'7FAAFF55';'7FD4FF2A'; ...
            '7FFFFF00';'7FFFD400';'7FFFAA00';'7FFF8000';'7FFF5500';'7FFF2A00';'7FFF0000';'50000000';'50000000']);
        
        fprintf(fida,'%s\r%s\r%s\r%s\r%s\r%s\r%s\r%s\r', ...
            line1,line2,line3,line4,line6,line7,line8,line9);
        
        %-----------------------------------
        % bin into value ranges (21 colours and 2 outside of range)
        %-----------------------------------
        
        deltaB=(maxBin-minBin)./20;
        f0=find(Tb<minBin);
        f1=find(Tb>=minBin & Tb<minBin+deltaB*1);
        f2=find(Tb>=minBin+deltaB*1 & Tb<minBin+deltaB*2);
        f3=find(Tb>=minBin+deltaB*2 & Tb<minBin+deltaB*3);
        f4=find(Tb>=minBin+deltaB*3 & Tb<minBin+deltaB*4);
        f5=find(Tb>=minBin+deltaB*4 & Tb<minBin+deltaB*5);
        f6=find(Tb>=minBin+deltaB*5 & Tb<minBin+deltaB*6);
        f7=find(Tb>=minBin+deltaB*6 & Tb<minBin+deltaB*7);
        f8=find(Tb>=minBin+deltaB*7 & Tb<minBin+deltaB*8);
        f9=find(Tb>=minBin+deltaB*8 & Tb<=minBin+deltaB*9);
        f10=find(Tb>=minBin+deltaB*9 & Tb<=minBin+deltaB*10);
        f11=find(Tb>=minBin+deltaB*10 & Tb<=minBin+deltaB*11);
        f12=find(Tb>=minBin+deltaB*11 & Tb<=minBin+deltaB*12);
        f13=find(Tb>=minBin+deltaB*12 & Tb<=minBin+deltaB*13);
        f14=find(Tb>=minBin+deltaB*13 & Tb<=minBin+deltaB*14);
        f15=find(Tb>=minBin+deltaB*14 & Tb<=minBin+deltaB*15);
        f16=find(Tb>=minBin+deltaB*15 & Tb<=minBin+deltaB*16);
        f17=find(Tb>=minBin+deltaB*16 & Tb<=minBin+deltaB*17);
        f18=find(Tb>=minBin+deltaB*17 & Tb<=minBin+deltaB*18);
        f19=find(Tb>=minBin+deltaB*18 & Tb<=minBin+deltaB*19);
        f20=find(Tb>=minBin+deltaB*19 & Tb<=minBin+deltaB*20);
        f21=find(Tb>maxBin);
        
        
        for j=1:22
            fprintf(fida,'%s %s %s %s %s\r','   <Placemark><name>TB Range: ',num2str((minBin+deltaB*(j-2))) ...
                ,'-',num2str((minBin+deltaB*(j-1))),'</name>');
            fprintf(fida,'%s\r','       <Style>');
            fprintf(fida,'%s\r','           <LineStyle>');
            fprintf(fida,'%s %s %s\r','               <color>',cls(j,:),'</color>');
            fprintf(fida,'%s\r','               <width> 0 </width>');
            fprintf(fida,'%s\r','           </LineStyle>');
            fprintf(fida,'%s\r','           <PolyStyle>');
            fprintf(fida,'%s %s %s\r','               <color>',cls(j,:),'</color>');
            fprintf(fida,'%s\r','           </PolyStyle>');
            fprintf(fida,'%s\r','       </Style>');
            fprintf(fida,'%s\r','     <MultiGeometry>');
            fbin = eval(genvarname(['f',num2str(j-1)]));
            if size(fbin ,1) > 1
                sizef = size(fbin,1);
            else
                sizef = size(fbin,2);
            end
            if ~isempty(fbin)
                k=fbin;
                for n=1:sizef
                    
                    fprintf(fida,'%s\r','       <Polygon>');
                    fprintf(fida,'%s\r','           <altitudeMode>relativetoground</altitudeMode>');
                    fprintf(fida,'%s\r','           <outerBoundaryIs>');
                    fprintf(fida,'%s\r','               <LinearRing>');
                    fprintf(fida,'%s\r','                   <coordinates>');
                    
                    for theta = 0:30:360
                        fprintf(fida,'%s%6f,%6f,%0f\r','                   ',Long(k(n))+footprint(k(n))*cosd(theta),Lat(k(n))+footprint(k(n))*sind(theta), 7000);
                    end
                    
                    
                    fprintf(fida,'%s\r','                   </coordinates>');
                    fprintf(fida,'%s\r','               </LinearRing>');
                    fprintf(fida,'%s\r','           </outerBoundaryIs>');
                    fprintf(fida,'%s\r','       </Polygon>');
                end
            end
            
            fprintf(fida,'%s\r','     </MultiGeometry>');
            fprintf(fida,'%s\r','   </Placemark>');
            
        end
        
        fprintf(fida,'%s\r',line_last);
        fclose(fida);
        
    end
end