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

% the azimuth scan angles were measured pos CCW from due North so we
% need to subtract them from 360 degrees to get them measured pos CW like
% heading and track angles which follow the aircraft campaign norms

az = 360 - azAll; 

az(az<0) = az(az<0) + 360;

 cycind = find(diff(az)>0)+1;
 tic
 offset = 0;
 for m = 1:length(az)
     if isempty(find(m == cycind, 1))
         az(m) = az(m) - offset;
     elseif ~isempty(find(m == cycind, 1))
         offset = offset + 360;
         az(m) = az(m) - offset;
     end
 end
 toc

azAll = az; 

save AzData timeMotAll azAll

quit

