%% SLAP Automated Report for the February 12, 2015 Flight

%% GPS Data Time History
% must run importOxTSPostprocessed.m first
clear
close all
% load GPS receiver data
oxfile = dir('OxTS*.mat');
load(oxfile(1).name)

% 16 timeond offset (converted to fraction of day)
% between the OxTS which is in GPS time and the rest of the
% data which is in UTC time
timeGeo = time + 16/86400;

% get more detail in timeGeo (raw data only specified to seconds but
% there are 100 data points with that second value in the time tag so the data
% freq is every 0.01 s. 
sgti = find(diff(timeGeo)>0.5/86400,1);
timeGeo2 = [ ((100 - sgti)/100/86400 + timeGeo(2)):0.01/86400:(timeGeo(end)+1/86400)]';

% add a nan in the first index of timeGeo because the OxTS data has a
% header line that must be accounted for
timeGeo = [nan; timeGeo2(1:length(timeGeo)-1)];

timeGeoInd = find(timeGeo >= timeGeo(2) & timeGeo < timeGeo(end) );
timeGeoSlice = timeGeo(timeGeoInd);
lat = PosLatdeg(timeGeoInd);
lon = PosLondeg(timeGeoInd);
alt = PosAltm(timeGeoInd);
roll = AngleRolldeg(timeGeoInd);
pitch = AnglePitchdeg(timeGeoInd);
hdg = AngleHeadingdeg(timeGeoInd);

% this section changes the trk data from being in cell format to vector
% format.
% trknew = AngleTrackdeg(timeGeoInd);
% trknew2 = nan(length(trknew),1);
% for c = 1:length(trknew)
%     trknew2(c) = str2double(trknew(c));
% end
% trk = trknew2;
% trk(trk<0) = trk(trk < 0) + 360;
% save trk_full trk

load trk_full

figure('units','normalized','outerposition',[0 0 1 1]) 
plot(lon, lat)
xlabel 'Longitude (deg)'
ylabel 'Latitude (deg)'

figure('units','normalized','outerposition',[0 0 1 1]) 
subplot(321)
plot(timeGeoSlice,alt), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Altitude (m)'
subplot(322)
plot(timeGeoSlice,roll), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Roll (deg)'
subplot(323)
plot(timeGeoSlice,hdg), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Heading (deg)'
subplot(324)
plot(timeGeoSlice,pitch), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Pitch (deg)'
subplot(325)
plot(timeGeoSlice,trk), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Track (deg)'
subplot(326)
yaw = hdg - trk;
plot(timeGeoSlice,yaw), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Yaw (deg)'

%% Motor Scan Angle Time History
if exist('AzData.mat')
    load AzData
else
    timeMotAll = [];
    azAll = [];
    filesENCODER = dir('ENC*');
    
    for m = 1:length(filesENCODER)
        
        [time,az_encoder,az_index] = import_encoder(filesENCODER(m).name);
        timeMot = time;
        az = mod(az_encoder - az_index,65536);
        
        % convert azimuth scan (az) angles that range from [0 65536] to be in range [0 360] deg
        az = az*360/65536;
        
        timeMotAll = [timeMotAll ; timeMot];
        azAll = [azAll ; az];

    end
    
    save AzData timeMotAll azAll
end

figure('units','normalized','outerposition',[0 0 1 1]) 

plot(timeMotAll, azAll), datetick
xlabel('Time (HH:MM UTC)');
ylabel 'Scan Angle (deg)'

%% Housekeeping Temperatures versus Time
if exist('daq_all.mat')
    load daq_all
else
    filesDaq = dir('ADAQTELEM*dlm');
    daq_all= [];
    for i = 1:length(filesDaq)
        temp = import_ADAQ_file(filesDaq(i).name);
        
        daq_all = [daq_all; temp];
    end
    save daq_all daq_all
end

timedaq = daq_all(:,1);
load hks_temp_names
figure('units','normalized','outerposition',[0 0 1 1])   
cnt = 0;
numfig = 1;
namecnt = 1;

limitflag = zeros(length(hks_names), 4);
for i = 2:57
    if cnt == 6
        figure('units','normalized','outerposition',[0 0 1 1]) 
        cnt = 0;
        numfig = numfig + 1;
    end
    cnt = cnt + 1;
    subplot(3,2,cnt)

    plot(daq_all(:,1), daq_all(:,i))

    ylabel('Temp (C)')
    xlabel('Time (HH:MM UTC)');
    datetick(gca)
    title(hks_names{namecnt}, 'fontsize', 12)
    namecnt = namecnt + 1;

end

%% radiometer data without noise
% the radiometer obesrvation data is named the same for both flights.  
% It is originally calculated in noise_remover.m for each 10 min
% data file then combined using the code shown below.
if exist('radiometer_wo_noise.mat')
    load radiometer_wo_noise
else
    filesRAD = dir('RAD*m2data.mat');
    [time_total, h_total, v_total]  = deal([]);
    for i = 10:length(filesRAD)-3
        
        filesRAD(i).name
        load(filesRAD(i).name)
        
        h_total = [h_total; h2new_v2];
        v_total = [v_total; v2new_v2];
        time_total = [time_total; time_file];
    end
    save radiometer_wo_noise time_total h_total v_total
end

figure('units','normalized','outerposition',[0 0 1 1]) 
subplot(211)
plot(time_total, h_total(1:length(time_total))), datetick
ylim([nanmean(h_total)-nanstd(h_total),nanmean(h_total)+nanstd(h_total)])
xlabel('Time (HH:MM UTC)');
ylabel 'Counts'
title 'H-Pol Antenna Observations'

subplot(212)
plot(time_total, v_total(1:length(time_total))), datetick
ylim([nanmean(v_total)-nanstd(v_total),nanmean(v_total)+nanstd(v_total)])
xlabel('Time (HH:MM UTC)');
ylabel 'Counts'
title 'V-Pol Antenna Observations'
% %% collect all reference data

% sample_period = 500e-6; %sec
% sample_period_days = sample_period/86400;
% 
% filesRAD = dir('RAD*FB_m2data.mat');
% [time_ref_total, h_ref_total, v_ref_total]  = deal([]);
% for i = 1:length(filesRAD)
%     
%     filesRAD(i).name
%     load(filesRAD(i).name)
%     
%     h_ref_total = [h_ref_total; h2ref];
%     v_ref_total = [v_ref_total; v2ref];
%     
%     time_file = [time, time + sample_period_days, time + sample_period_days*2, time + sample_period_days*3];
%     time_file =  reshape(time_file', numel(time_file), 1);
%     
%     time_ref_total = [time_ref_total; time_file];
% end
% 
% ref_nonnan = find(~isnan(h_ref_total));
% 
% figure('units','normalized','outerposition',[0 0 1 1]) 
% subplot(211)
% plot(time_ref_total(ref_nonnan), h_ref_total(ref_nonnan),'.'), datetick
% ylim([nanmean(h_ref_total)-nanstd(h_ref_total),nanmean(h_ref_total)+nanstd(h_ref_total)])
% xlabel('Time (HH:MM UTC)');
% ylabel 'Counts'
% title 'H-Pol Reference Observations'
% 
% subplot(212)
% plot(time_ref_total(ref_nonnan), v_ref_total(ref_nonnan),'.'), datetick
% ylim([nanmean(v_ref_total)-nanstd(v_ref_total),nanmean(v_ref_total)+nanstd(v_ref_total)])
% xlabel('Time (HH:MM UTC)');
% ylabel 'Counts'
% title 'V-Pol Reference Observations'


