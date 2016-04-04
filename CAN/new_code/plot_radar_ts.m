%%
clear

%% Initialize key options
flagHalf = 1;
% flag to process data for fore or aft half-scan
fore = 1;

% flag to plot h or v pol data on Google Earth
GoogleEarth  = 1;

%% roll flag and limits
flagRoll = 1;
% Aircraft roll angle limits if removing data during a turn
rollmin = -5;
rollmax = 5;

% RDRRETURN_20151030T181002.mat
files = dir('6deg_RDRRETURN_*.mat');
files_to_process = 1:length(files);
%%
% beam width information by polarization from SLAP documentation
beam_width_h_pol = 20.5;
beam_width_v_pol = 19.5;

[time_total, az_total, lat_total, lon_total, hh_total, vv_total, hv_total, ...
 vh_total, alt_total, trk_total, roll_total, hdg_total] = deal([]);


%%
for k = files_to_process
    files(k).name
    % load time, radar_h, and radar_v data
    load(files(k).name)
    
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

end   % file loop 

figure
subplot(3, 1, 1)
plot(time_total, vv_total, 'b')
title('Time series of radar sigma0 VV (dB)')
axis([-inf, inf, -20, 20])
datetick('x', 15)

subplot(3, 1, 2)
plot(time_total, hh_total, 'g')
axis([-inf, inf, -20, 20])
title('Time series of radar sgima0 HH (dB)')
datetick('x', 15)

subplot(3, 1, 3)
plot(time_total, roll_total)
title('Roll angle (deg) radar')
datetick('x', 15)

print('radar_vv_hh_ts.png', '-dpng');
