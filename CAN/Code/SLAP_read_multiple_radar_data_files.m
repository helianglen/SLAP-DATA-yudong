% function SLAP_read_multiple_radar_data_files(val, loopback)

val = 1:3;

% if nargin < 2
    files = dir('RDRRET*.slapbin');
% else 
%     files = dir('RDRLOOP*.slapbin');
% end

for i = val
    SLAP_read_raw_radar_data(files(i).name);
end