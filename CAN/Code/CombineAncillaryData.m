%% import then combine all daq data
filesDaq = dir('ADAQTELEM*dlm');
daq_all= [];
for i = 1:length(filesDaq)
    temp = import_ADAQ_file(filesDaq(i).name);
    
    daq_all = [daq_all; temp];
end
    save daq_all daq_all
    
%% combine all daq data if imported data already exists (use with iPHEx campaign data)
filesDaq = dir('DAQTELEM*2.mat');
daq_all= [];
for i = 1:length(filesDaq)
    load(filesDaq(i).name)
    daq_all = [daq_all; temp(:,1:28)];
end
save daq_all daq_all
%% combine all encoder scan angle data
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

%% combine all mot data
filesMot = dir('MOT*slapdlm');

[timeMotAll, resolver_all, motor_current_all]  = deal([]);
for i = 11:21%1:length(filesMot)
    filesMot(i).name
    [time,resolver, motor_current] = import_motion(filesMot(i).name);
    timeMotAll = [timeMotAll; datenum(time)];
    resolver_all = [resolver_all; resolver];
    motor_current_all = [motor_current_all; motor_current];
end
% convert az angles that range from [0 65535] to be in range [0 360]
azAll = resolver_all*360/65536;

save RESOLVER_AzData timeMotAll azAll

% figure; plot(time_mot, motor_current_all, '.')
% datetick(gca)
% title 'May 21 AM'
% xlabel 'Time'
% ylabel 'Motor Current'
% saveas(gcf, 'may 21 AM motor current.jpg')
    