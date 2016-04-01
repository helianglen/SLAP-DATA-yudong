%function SLAP_FBM2_geolocate(data_path)
clc
clear all
close all

tic

name = '20151106T124910 SLAPex CAN Morning Flight';

data_path = '.';  
mean_elev = 238; %m

ingest_script_path = pwd;
GE_script_path = fullfile(pwd, 'SLAP_GE_plot/');

%% Read in radar data and split the poles

fprintf('\n');
disp('Ingesting Radar Data')

cd(data_path);
%data_files = dir('RDRRETURN*.mat');
data_files = dir('RDRRETURN_20151030T185000.mat');
cd(ingest_script_path);

[radar_hh, radar_hv, radar_vh, radar_vv, radar_t] = deal([]);

for i = 1:length(data_files)
    
    disp(['Reading file: ' data_files(i).name '...'])
    
    load(fullfile(data_path, data_files(i).name))
    
    [radar_t_temp, radar_hh_temp, radar_hv_temp, radar_vh_temp, radar_vv_temp] = SLAP_radar_split_poles(time, radar_h, radar_v, pcm);
    
    radar_hh = [radar_hh; radar_hh_temp];
    radar_hv = [radar_hv; radar_hv_temp];
    radar_vh = [radar_vh; radar_vh_temp];
    radar_vv = [radar_vv; radar_vv_temp];
    
    radar_t = [radar_t; radar_t_temp];

end

keyboard

%% Read in position/attitude data and interpolate

fprintf('\n');
disp('Ingesting AHRS Data')

cd(data_path);
ahrs_files = dir('AHRSTELEM*.slapdlm');
%cd(ingest_script_path);

[ahrs_t, ahrs_lat, ahrs_lon, ahrs_alt, ahrs_roll, ahrs_pitch, ahrs_heading, ahrs_gps_t, ahrs_gps_t_valid] = deal([]);

for i = 1:length(ahrs_files)
    
    disp(['Reading file: ' ahrs_files(i).name '...'])
    
    [ahrs_t_temp, ahrs_lat_temp, ahrs_lon_temp, ahrs_alt_temp, ~, ~, ~, ahrs_roll_temp, ahrs_pitch_temp, ahrs_heading_temp, ahrs_gps_t_temp, ahrs_gps_t_valid_temp, ~, ~, ~, ~,~, ~, ~,~, ~, ~, ~] = import_ahrs(fullfile(data_path, ahrs_files(i).name));
    
    ahrs_t = [ahrs_t; ahrs_t_temp];
    ahrs_lat = [ahrs_lat; ahrs_lat_temp];
    ahrs_lon = [ahrs_lon; ahrs_lon_temp];
    ahrs_alt = [ahrs_alt; ahrs_alt_temp];
    ahrs_roll = [ahrs_roll; ahrs_roll_temp];
    ahrs_pitch = [ahrs_pitch; ahrs_pitch_temp];
    ahrs_heading = [ahrs_heading; ahrs_heading_temp];
    ahrs_gps_t = [ahrs_gps_t; ahrs_gps_t_temp];
    ahrs_gps_t_valid = [ahrs_gps_t_valid; ahrs_gps_t_valid_temp];
    
end

% For some reason the AHRS timestamps aren't monotonically increasing.
% Until we can get the post-processed OxTS data, sort the data by the
% timestamps.

[~, ahrs_sort_i] = sort(ahrs_gps_t);

% Remove invalid AHRS data
ahrs_sort_i = ahrs_sort_i(ahrs_gps_t_valid == 1);

% Find repeat points and remove
ahrs_sort_i = ahrs_sort_i([diff(ahrs_gps_t(ahrs_sort_i)) ~= 0; true]);

ahrs_lat_interp = interp1(ahrs_gps_t(ahrs_sort_i), ahrs_lat(ahrs_sort_i), radar_t);
ahrs_lon_interp = interp1(ahrs_gps_t(ahrs_sort_i), ahrs_lon(ahrs_sort_i), radar_t);
ahrs_alt_interp = interp1(ahrs_gps_t(ahrs_sort_i), ahrs_alt(ahrs_sort_i), radar_t);
ahrs_roll_interp = interp1(ahrs_gps_t(ahrs_sort_i), ahrs_roll(ahrs_sort_i), radar_t);
ahrs_pitch_interp = interp1(ahrs_gps_t(ahrs_sort_i), ahrs_pitch(ahrs_sort_i), radar_t);
ahrs_heading_interp = interp1(ahrs_gps_t(ahrs_sort_i), ahrs_heading(ahrs_sort_i), radar_t);

ahrs_agl = ahrs_alt_interp - mean_elev;

%% Read in azimuth data and interpolate

fprintf('\n');
disp('Ingesting Encoder Data')

cd(data_path);
enc_files = dir('ENCODERTELEM*.slapdlm');
%cd(ingest_script_path);

[enc_t, enc_az, enc_az_index] = deal([]);

for i = 1:length(enc_files)
    
    disp(['Reading file: ' enc_files(i).name '...'])
    
    [enc_t_temp, enc_az_temp, enc_az_index_temp, ~, ~] = import_encoder(fullfile(data_path, enc_files(i).name));
    
    enc_t = [enc_t; enc_t_temp];
    enc_az = [enc_az; enc_az_temp];
    enc_az_index = [enc_az_index; enc_az_index_temp];
    
end

% Sort encoder data by timestamp and remove duplicates
[~, enc_sort_i] = sort(enc_t);
enc_sort_i = enc_sort_i([diff(enc_t(enc_sort_i)) ~= 0; true]);

enc_az_interp = interp1(enc_t(enc_sort_i), enc_az(enc_sort_i), radar_t);
enc_az_index_interp = interp1(enc_t(enc_sort_i), enc_az_index(enc_sort_i), radar_t, 'nearest');

enc_counts = mod(enc_az_interp - enc_az_index_interp, 65536);
enc_deg = enc_counts/65536*360;

%% Compute geolocation for data points

fprintf('\n');
disp('Computing geolocation...');

% Constants
phi = 40; %deg - SLAP nadir angle
r_earth = 6371e3; %m - Radius of the Earth

% Radius of conical scan
scan_radius = tand(phi) * ahrs_agl;

% Slant range
slant = ahrs_agl/cosd(phi);

% Calculate beam heading
beam_heading = ahrs_heading_interp - enc_deg;

lat_offset_m = scan_radius .* cosd(beam_heading); %m
lon_offset_m = scan_radius .* sind(beam_heading); %m

lat_offset_deg = lat_offset_m * 360/(2*pi*r_earth); %deg
lon_offset_deg = lon_offset_m * 360/(2*pi*r_earth); %deg

data_lat = ahrs_lat_interp + lat_offset_deg;
data_lon = ahrs_lon_interp + lon_offset_deg;

beamwidth = 12; %deg

% new equation to determine footprint size
footprintsm = (slant*tand(beamwidth/2))/cosd(phi);

% convert footprint radius in meters to degrees longitude
footprints = footprintsm / (pi*r_earth) * 180;

%% Integrate radar power

fprintf('\n');
disp('Integrating range bins...');

numbins = size(radar_h, 2);

%first radar return bin is at 4 us after transmit pulse
rangegate = ((0:numbins-1) + 4) * 1e-6; %sec

c = 299792458; %m/s

range_multiplier = repmat(rangegate * c, size(radar_hh, 1), 1);

radar_hh_sum = sum(radar_hh .* (range_multiplier.^4), 2);
radar_hv_sum = sum(radar_hv .* (range_multiplier.^4), 2);
radar_vh_sum = sum(radar_vh .* (range_multiplier.^4), 2);
radar_vv_sum = sum(radar_vv .* (range_multiplier.^4), 2);

'check range' 
keyboard

%% Average data

fprintf('\n');
disp('Averaging data...');

% For a cross-scan beam-width of 12 deg, average all data in a 6 degree
% swath for spatial nyquist. This works out to 132 consecutive samples at
% 15 RPM with a 500 microsecond PRI
%
% (132 * 500 microseconds * 15 (360 degrees / min) in degrees = 5.94 deg)

avg_samples = 132;

num_averages = floor(length(radar_t)/avg_samples);

[data_lat_avg, data_lon_avg, radar_hh_avg, radar_hv_avg, radar_vh_avg, radar_vv_avg, ahrs_roll_avg, footprints_avg] = deal(nan(num_averages, 1));

for i = 1:num_averages
    avg_start = (i-1)*avg_samples+1;
    avg_end = i*avg_samples;
    
    data_lat_avg(i) = nanmean(data_lat(avg_start:avg_end));
    data_lon_avg(i) = nanmean(data_lon(avg_start:avg_end));
    
    radar_hh_avg(i) = nanmean(radar_hh_sum(avg_start:avg_end));
    radar_hv_avg(i) = nanmean(radar_hv_sum(avg_start:avg_end));
    radar_vh_avg(i) = nanmean(radar_vh_sum(avg_start:avg_end));
    radar_vv_avg(i) = nanmean(radar_vv_sum(avg_start:avg_end));
    
    ahrs_roll_avg(i) = nanmean(ahrs_roll_interp(avg_start:avg_end));
    footprints_avg(i) = nanmean(footprints(avg_start:avg_end));
    
end

%% Convert to dB

%remove zeroes/negatives before conversion to dB
radar_hh_avg(radar_hh_avg <= 0) = nan;
radar_hv_avg(radar_hv_avg <= 0) = nan;
radar_vh_avg(radar_vh_avg <= 0) = nan;
radar_vv_avg(radar_vv_avg <= 0) = nan;

radar_hh_dB = 10*log10(radar_hh_avg);
radar_hv_dB = 10*log10(radar_hv_avg);
radar_vh_dB = 10*log10(radar_vh_avg);
radar_vv_dB = 10*log10(radar_vv_avg);

%% Filter data

fprintf('\n');
disp('Filtering data...');

max_roll = 5; %deg

% Filling in NaN for the latitude will prevent the point from being plotted
data_lat_avg(abs(ahrs_roll_avg) > 5) = nan;

%% Compute color bins

fprintf('\n');
disp('Computing bins...');

cd(GE_script_path);

radar_hh_mean = nanmean(radar_hh_dB);
%radar_hh_std = std(radar_hh_dB, 'omitnan');
radar_vv_mean = nanmean(radar_vv_dB);
%radar_vv_std = std(radar_vv_dB, 'omitnan');

%MATLAB 2013 compatibility
radar_hh_nonan = radar_hh_dB(~isnan(radar_hh_dB));
radar_vv_nonan = radar_vv_dB(~isnan(radar_vv_dB));

radar_hh_std = std(radar_hh_nonan);
radar_vv_std = std(radar_vv_nonan);

%Top of bottom bin is 2 std below mean, bottom of top bin is 2 std above
%mean
radar_hh_min = radar_hh_mean - 2*radar_hh_std;
radar_hh_max = radar_hh_mean + 2*radar_hh_std;

radar_vv_min = radar_vv_mean - 2*radar_vv_std;
radar_vv_max = radar_vv_mean + 2*radar_vv_std;

[ radar_hh_lat_bins, radar_hh_lon_bins, radar_hh_counts_bins, radar_hh_footprints_bins, radar_hh_beam_heading_bins ] = GE_gen_color_bins( data_lat_avg, data_lon_avg, footprints_avg, radar_hh_dB, radar_hh_min, radar_hh_max, beam_heading, 22 );
[ radar_vv_lat_bins, radar_vv_lon_bins, radar_vv_counts_bins, radar_vv_footprints_bins, radar_vv_beam_heading_bins ] = GE_gen_color_bins( data_lat_avg, data_lon_avg, footprints_avg, radar_vh_dB, radar_vv_min, radar_vv_max, beam_heading, 22 );

%% Generate KML file

fprintf('\n');
disp('Generating KML files...');

GE_generate_kml( data_path, [name '_radar_hh.kml'], [name ' Radar HH'], radar_hh_lat_bins, radar_hh_lon_bins, radar_hh_footprints_bins, radar_hh_footprints_bins, radar_hh_counts_bins, radar_hh_beam_heading_bins );
GE_generate_kml( data_path, [name '_radar_vv.kml'], [name ' Radar VV'], radar_vv_lat_bins, radar_vv_lon_bins, radar_vv_footprints_bins, radar_vv_footprints_bins, radar_vv_counts_bins, radar_vv_beam_heading_bins );

%cd(ingest_script_path);

elapsed = toc;
disp(['Time elapsed: ' num2str(elapsed) ' s.']);
