%% combine all encoder scan angle data
timeMotAll = [];
azAll = [];

% tempariy variables 
[enc_t, enc_az, enc_az_index] = deal([]);

filesENCODER = dir('ENCODERTELEM*.slapdlm');

for m = 1:length(filesENCODER)

    [time,az_encoder,az_index] = import_encoder(filesENCODER(m).name);
   
    enc_t = [enc_t; time]; 
    enc_az = [enc_az; az_encoder];
    enc_az_index = [enc_az_index; az_index]; 

end 

% Albert's code: sort encoder data by timestamp and remove duplicates
[~, enc_sort_i] = sort(enc_t);
enc_sort_i = enc_sort_i([diff(enc_t(enc_sort_i)) ~= 0; true]);

timeMotAll = enc_t(enc_sort_i);  
azAll = mod(enc_az(enc_sort_i) -  enc_az_index(enc_sort_i, 65535);
azAll = azAll*360/65536; 

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

