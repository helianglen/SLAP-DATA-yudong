clear
files = dir('RDRRET*6deg.mat');

[time_total, az_total, lat_total, lon_total, hh_total, vv_total, hv_total, vh_total, alt_total, ...
    trk_total, roll_total, hdg_total] = deal([]);
for i = 1:length(files)
    
    load(files(i).name)
    time_total = [time_total; timeavg];
    hh_total = [hh_total; mean_hh];
    vv_total = [vv_total; mean_vv];
    hv_total = [hv_total; mean_hv];
    vh_total = [vh_total; mean_vh];
    az_total = [az_total; azavg];
    lat_total = [lat_total; latavg];
    lon_total = [lon_total; lonavg];
    alt_total = [alt_total; altavg];
    trk_total = [trk_total; trkavg];
    roll_total = [roll_total; rollavg];
    hdg_total = [hdg_total; hdgavg];
end

% set flag to determine if SLAP is spinning or not based on the
% difference in consecutive azimuth scan angles being less than 4 deg
scan_motor_diff = abs(diff(az_total));
scan_motor_spinning = scan_motor_diff > 4;
scan_motor_spinning = [scan_motor_spinning; logical(0)];

% determine if data point is in the fore or aft half scan
% (N/A to stare data)
flagHS = az_total > 270 | az_total < 90;

% combine all data into one matrix
alldata = [ time_total, Lat_total, Lon_total, alt_total, tbh_total, tbv_total,  az_total, roll_total, trk_total, hdg_total, scan_motor_spinning, flagHS];

% remove data lines with all NaN values
indnan = isnan(Lat_total);
alldata = alldata(~indnan,:);

% save as csv file
csvwrite([ 'SLAP May 21 PM Radar Observations.csv'] , alldata)


% figure; plot(time_total, h_total, '.')
% hold on; plot(time_total, v_total, '.')
% datetick