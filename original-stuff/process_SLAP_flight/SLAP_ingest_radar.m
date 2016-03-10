function SLAP_ingest_radar(filename, numpackets)

tic

%Check endianness for typecasting. If system is little-endian, swap bytes
%of big-endian FPGA data
[~, ~, endianness] = computer;
swap_endianness = endianness == 'L';


disp(filename);

fh = fopen(filename);
filedata = fread(fh);
fclose(fh);

numbytes = length(filedata);

if numbytes < 4
    disp('File too small.');
    return
end

bytesperpacket = SLAP_UINT32(filedata(1:4), swap_endianness); %read the first size field to get the size of all bytes, depends on the radar receive window, but should be around 174 for the return and 34-70 for the loopback
packetid_start = 1;
opmode_start = 2;
time_start = 3;
pcm_start = 19;
data_start = 23;
data_size = bytesperpacket - 22;

if numbytes < bytesperpacket
    return
end

first_time_bytes = filedata(7:22);
first_time = SLAP_IRIG_ns_since_TAI_v2(first_time_bytes, swap_endianness);
try
    disp(['First packet at ' datestr(first_time)]);
catch err
    if strcmp(err.identifier, 'MATLAB:datestr:ConvertDateNumber')
        disp(['Warning: First timestamp invalid - check file integrity (timestamp: ' num2str(first_time_bytes', '%02x') ')']);
    end
end

if nargin < 2
    numpackets = floor(numbytes/(bytesperpacket+4));
end

disp(['Processing ' num2str(numpackets) ' packets...']);

[radar_h,   radar_v] = deal(nan(numpackets, data_size/4));

[mode, time, pcm, pulses, cycles] = deal(nan(numpackets,1));
time_bytes = nan(numpackets, 16);

filedata = reshape(filedata(1:(bytesperpacket+4)*numpackets), bytesperpacket+4, numpackets)';

parfor i=1:numpackets
    
    packet = filedata(i,:);
    packet = packet(5:end);
    
    mode(i) = packet(opmode_start);
    
    time_bytes(i,:) = packet(time_start:time_start+16-1);
    [time(i), pulses(i), cycles(i)] = SLAP_IRIG_ns_since_TAI_v2(time_bytes(i,:), swap_endianness);
    
    pcm(i) = SLAP_UINT32(packet(pcm_start:pcm_start+4-1), swap_endianness);
    
    rawdat = packet(data_start:data_start+data_size-1);
    
    %% Radar data
    
    [h, v] = deal(nan(1, data_size/4));
    
    v_msb = rawdat(1:4:end);
    v_lsb = rawdat(2:4:end);
    h_msb = rawdat(3:4:end);
    h_lsb = rawdat(4:4:end);
    
    
    for j=1:data_size/4
        h(j) = SLAP_UINT16([h_msb(j) h_lsb(j)], swap_endianness);
        v(j) = SLAP_UINT16([v_msb(j) v_lsb(j)], swap_endianness);
    end
    
    radar_h(i, :) = h;
    radar_v(i, :) = v;
    
%     if sum(i == 1:floor(numpackets/10):numpackets)
%         disp([num2str(i), 'out of ', num2str(numpackets), 'completed!'])
%     end
end

time = SLAP_fix_timestamps(time, pulses, cycles, 500e-6);

elapsedtime = toc;

disp(['Execution Time: ' num2str(elapsedtime) ' s']);

[pathstr,name,ext] = fileparts(filename);
filename_wo_ext = fullfile(pathstr,name);
save(filename_wo_ext, 'time', 'pcm', 'mode', 'radar_h', 'radar_v')

