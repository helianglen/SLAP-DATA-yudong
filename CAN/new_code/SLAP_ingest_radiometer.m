function SLAP_ingest_radiometer(filename, varargin)
%SLAP_L1B_32 Parse SLAP radiometer output .SLAPBIN file and save as .MAT
%   SLAP_L1B_32(FILENAME) parses fullband data from a SLAP radiometer
%   output binary file into a MATLAB formatted binary file (.MAT).
%
%   SLAP_L1B_32(FILENAME, OPTIONS) parses fullband and/or subband data from
%   a SLAP radiometer output binary file into a MATLAB formatted binary
%   file (.MAT).
%
%   Inputs:
%
%   FILENAME: The path and filename of the .SLAPBIN file containing 32-bit
%   SLAP radiometer data.
%
%   OPTIONS: Specify whether to process fullband data, subband data, or
%   both. If neither is specified, only fullband data will be parsed. Use
%   one of the following combinations:
%
%       'fb'                        Only parse fullband data (default).
%       'sb'                        Only parse subband data.
%       'fb', 'sb'                  Parse both fullband and subband data.

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

bytesperpacket = SLAP_UINT32(filedata(1:4), swap_endianness); %read the first size field to get the size of all bytes, should be 1142
packetid_start = 1;
opmode_start = 2;
time_start = 3;
pcm_start = 19;
data_start = 23;
data_size = 1120;

if numbytes < bytesperpacket
    return
end

if bytesperpacket ~= 1142
    disp('Packet size invalid for 32-bit SLAP radiometer format. Please try a different file.');
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

numpackets = floor(numbytes/(bytesperpacket+4));

process_FB = 0;
process_SB = 0;

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i}, 'FB')
            process_FB = 1;
        elseif strcmpi(varargin{i}, 'SB')
            process_SB = 1;
        end
    end
end

if ~(process_FB || process_SB)
    process_FB = 1;
end


disp(['Processing ' num2str(numpackets) ' packets...']);
if process_FB
    [m1_ant,     m2_ant,      m3_ant,          m4_ant,      ...
        m1_ant_nd,  m2_ant_nd,   m3_ant_nd,       m4_ant_nd,   ...
        m1_ref,     m2_ref,      m3_ref,          m4_ref,      ...
        m1_ref_nd,  m2_ref_nd,   m3_ref_nd,       m4_ref_nd] = deal(nan(numpackets, 4, 4));
    
    [t3_ant,     t3_ant_nd,   t3_ref,          t3_ref_nd,   ...
        t4_ant,     t4_ant_nd,   t4_ref,          t4_ref_nd] = deal(nan(numpackets, 4));
end
if process_SB
    [m1_16_ant,  m1_16_ref,   m1_16_ant_nd,    m1_16_ref_nd,...
        m2_16_ant,  m2_16_ref,   m2_16_ant_nd,    m2_16_ref_nd,...
        m3_16_ant,  m3_16_ref,   m3_16_ant_nd,    m3_16_ref_nd,...
        m4_16_ant,  m4_16_ref,   m4_16_ant_nd,    m4_16_ref_nd] = deal(nan(numpackets, 4, 16));
    
    [t3_16_ant,  t3_16_ref,   t3_16_ant_nd,    t3_16_ref_nd,...
        t4_16_ant,  t4_16_ref,   t4_16_ant_nd,    t4_16_ref_nd] = deal(nan(numpackets, 16));
end

[mode, time, pulses, cycles, pcm] = deal(nan(numpackets,1));
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
    %% Subband data
    if process_SB
        
        sb_bytes = rawdat(1:768);
        sb_bytes = reshape(sb_bytes, 48, 16)';
        
        [subHM1I, subVM1I, subHM1Q, subVM1Q, ...
            subHM2I, subVM2I, subHM2Q, subVM2Q, ...
            subHM3I, subVM3I, subHM3Q, subVM3Q, ...
            subHM4I, subVM4I, subHM4Q, subVM4Q] = deal(nan(1, 16));
        
        for j = 1:16
            subHM1I(j) = SLAP_INT16(sb_bytes(j, 1:2), swap_endianness);
            subHM2I(j) = SLAP_UINT32(sb_bytes(j, 5:8), swap_endianness);
            subHM3I(j) = SLAP_INT16(sb_bytes(j, 13:14), swap_endianness);
            subHM4I(j) = SLAP_UINT32(sb_bytes(j, 17:20), swap_endianness);
            
            subVM1I(j) = SLAP_INT16(sb_bytes(j, 3:4), swap_endianness);
            subVM2I(j) = SLAP_UINT32(sb_bytes(j, 9:12), swap_endianness);
            subVM3I(j) = SLAP_INT16(sb_bytes(j, 15:16), swap_endianness);
            subVM4I(j) = SLAP_UINT32(sb_bytes(j, 21:24), swap_endianness);
            
            subHM1Q(j) = SLAP_INT16(sb_bytes(j, 25:26), swap_endianness);
            subHM2Q(j) = SLAP_UINT32(sb_bytes(j, 29:32), swap_endianness);
            subHM3Q(j) = SLAP_INT16(sb_bytes(j, 37:38), swap_endianness);
            subHM4Q(j) = SLAP_UINT32(sb_bytes(j, 41:44), swap_endianness);
            
            subVM1Q(j) = SLAP_INT16(sb_bytes(j, 27:28), swap_endianness);
            subVM2Q(j) = SLAP_UINT32(sb_bytes(j, 33:36), swap_endianness);
            subVM3Q(j) = SLAP_INT16(sb_bytes(j, 39:40), swap_endianness);
            subVM4Q(j) = SLAP_UINT32(sb_bytes(j, 45:48), swap_endianness);
        end
    end
    
    %% Fullband Data
    if process_FB
        fb_bytes = rawdat(769:960);
        fb_bytes = reshape(fb_bytes, 48, 4)';
        
        [HM1I, VM1I, HM1Q, VM1Q, ...
            HM2I, VM2I, HM2Q, VM2Q, ...
            HM3I, VM3I, HM3Q, VM3Q, ...
            HM4I, VM4I, HM4Q, VM4Q] = deal(nan(1, 4));
        
        for j = 1:4
            HM1I(j) = SLAP_INT16(fb_bytes(j, 1:2), swap_endianness);
            HM2I(j) = SLAP_UINT32(fb_bytes(j, 9:12), swap_endianness);
            HM3I(j) = SLAP_INT16(fb_bytes(j, 25:26), swap_endianness);
            HM4I(j) = SLAP_UINT32(fb_bytes(j, 33:36), swap_endianness);
            
            VM1I(j) = SLAP_INT16(fb_bytes(j, 3:4), swap_endianness);
            VM2I(j) = SLAP_UINT32(fb_bytes(j, 13:16), swap_endianness);
            VM3I(j) = SLAP_INT16(fb_bytes(j, 27:28), swap_endianness);
            VM4I(j) = SLAP_UINT32(fb_bytes(j, 37:40), swap_endianness);
            
            HM1Q(j) = SLAP_INT16(fb_bytes(j, 5:6), swap_endianness);
            HM2Q(j) = SLAP_UINT32(fb_bytes(j, 17:20), swap_endianness);
            HM3Q(j) = SLAP_INT16(fb_bytes(j, 29:30), swap_endianness);
            HM4Q(j) = SLAP_UINT32(fb_bytes(j, 41:44), swap_endianness);
            
            VM1Q(j) = SLAP_INT16(fb_bytes(j, 7:8), swap_endianness);
            VM2Q(j) = SLAP_UINT32(fb_bytes(j, 21:24), swap_endianness);
            VM3Q(j) = SLAP_INT16(fb_bytes(j, 31:32), swap_endianness);
            VM4Q(j) = SLAP_UINT32(fb_bytes(j, 45:48), swap_endianness);
        end
    end
    
    %% Subband Crosscorrelation
    if process_SB
        sbxc_bytes = rawdat(961:1088);
        
        [subRHVI, subRHVQ] = deal(nan(1, 16));
        
        for j = 1:16
            subRHVI(j) = SLAP_UINT32(sbxc_bytes((j-1)*8+1:j*8-4), swap_endianness);
            subRHVQ(j) = SLAP_UINT32(sbxc_bytes((j-1)*8+5:j*8), swap_endianness);
        end
    end
    
    %% Fullband Crosscorrelation
    if process_FB
        fbxc_bytes = rawdat(1089:end);
        
        [RHVI, RHVQ] = deal(nan(1, 4));
        
        for j = 1:4
            RHVI(j) = SLAP_UINT32(fbxc_bytes((j-1)*8+1:j*8-4), swap_endianness);
            RHVQ(j) = SLAP_UINT32(fbxc_bytes((j-1)*8+5:j*8), swap_endianness);
        end
    end
    
    %% Sort based on PCM bitfield
    
    if bitget(pcm(i), 7) && bitget(pcm(i), 8) % ref, no INS
        
        %fullband indexing: mx_(packet, polarization/moment, PRI)
        if process_FB
            m1_ref(i,:,:)       = [HM1I; HM1Q; VM1I; VM1Q];
            m2_ref(i,:,:)       = [HM2I; HM2Q; VM2I; VM2Q];
            m3_ref(i,:,:)       = [HM3I; HM3Q; VM3I; VM3Q];
            m4_ref(i,:,:)       = [HM4I; HM4Q; VM4I; VM4Q];
            t3_ref(i,:)         = RHVI;
            t4_ref(i,:)         = RHVQ;
        end
        
        %subband indexing: mx_16_(packet, polarization/moment, subband)
        if process_SB
            m1_16_ref(i,:,:)    = [subHM1I; subHM1Q; subVM1I; subVM1Q];
            m2_16_ref(i,:,:)	= [subHM2I; subHM2Q; subVM2I; subVM2Q];
            m3_16_ref(i,:,:)	= [subHM3I; subHM3Q; subVM3I; subVM3Q];
            m4_16_ref(i,:,:)	= [subHM4I; subHM4Q; subVM4I; subVM4Q];
            t3_16_ref(i,:,:)	= subRHVI;
            t4_16_ref(i,:,:)	= subRHVQ;
        end
        
        
    elseif bitget(pcm(i), 7) % ref, INS
        if process_FB
            m1_ref_nd(i,:,:)	= [HM1I; HM1Q; VM1I; VM1Q];
            m2_ref_nd(i,:,:)	= [HM2I; HM2Q; VM2I; VM2Q];
            m3_ref_nd(i,:,:)	= [HM3I; HM3Q; VM3I; VM3Q];
            m4_ref_nd(i,:,:)    = [HM4I; HM4Q; VM4I; VM4Q];
            t3_ref_nd(i,:)      = RHVI;
            t4_ref_nd(i,:)      = RHVQ;
        end
        if process_SB
            m1_16_ref_nd(i,:,:) = [subHM1I; subHM1Q; subVM1I; subVM1Q];
            m2_16_ref_nd(i,:,:)	= [subHM2I; subHM2Q; subVM2I; subVM2Q];
            m3_16_ref_nd(i,:,:)	= [subHM3I; subHM3Q; subVM3I; subVM3Q];
            m4_16_ref_nd(i,:,:)	= [subHM4I; subHM4Q; subVM4I; subVM4Q];
            t3_16_ref_nd(i,:,:)	= subRHVI;
            t4_16_ref_nd(i,:,:)	= subRHVQ;
        end
        
    elseif bitget(pcm(i), 8) % antenna, no INS
        if process_FB
            m1_ant(i,:,:)       = [HM1I; HM1Q; VM1I; VM1Q];
            m2_ant(i,:,:)       = [HM2I; HM2Q; VM2I; VM2Q];
            m3_ant(i,:,:)       = [HM3I; HM3Q; VM3I; VM3Q];
            m4_ant(i,:,:)       = [HM4I; HM4Q; VM4I; VM4Q];
            t3_ant(i,:)         = RHVI;
            t4_ant(i,:)         = RHVQ;
        end
        if process_SB
            m1_16_ant(i,:,:)    = [subHM1I; subHM1Q; subVM1I; subVM1Q];
            m2_16_ant(i,:,:)	= [subHM2I; subHM2Q; subVM2I; subVM2Q];
            m3_16_ant(i,:,:)	= [subHM3I; subHM3Q; subVM3I; subVM3Q];
            m4_16_ant(i,:,:)	= [subHM4I; subHM4Q; subVM4I; subVM4Q];
            t3_16_ant(i,:,:)	= subRHVI;
            t4_16_ant(i,:,:)	= subRHVQ;
        end
    else % antenna, INS
        if process_FB
            m1_ant_nd(i,:,:)	= [HM1I; HM1Q; VM1I; VM1Q];
            m2_ant_nd(i,:,:)	= [HM2I; HM2Q; VM2I; VM2Q];
            m3_ant_nd(i,:,:)	= [HM3I; HM3Q; VM3I; VM3Q];
            m4_ant_nd(i,:,:)    = [HM4I; HM4Q; VM4I; VM4Q];
            t3_ant_nd(i,:)      = RHVI;
            t4_ant_nd(i,:)      = RHVQ;
        end
        if process_SB
            m1_16_ant_nd(i,:,:) = [subHM1I; subHM1Q; subVM1I; subVM1Q];
            m2_16_ant_nd(i,:,:)	= [subHM2I; subHM2Q; subVM2I; subVM2Q];
            m3_16_ant_nd(i,:,:)	= [subHM3I; subHM3Q; subVM3I; subVM3Q];
            m4_16_ant_nd(i,:,:)	= [subHM4I; subHM4Q; subVM4I; subVM4Q];
            t3_16_ant_nd(i,:,:)	= subRHVI;
            t4_16_ant_nd(i,:,:)	= subRHVQ;
        end
    end

end

time = SLAP_fix_timestamps(time, pulses, cycles, 2000e-6); %fix timestamps

if process_FB
    FullMomGroup = struct('m1_ant', m1_ant, 'm1_ant_nd', m1_ant_nd, 'm1_ref', m1_ref, 'm1_ref_nd', m1_ref_nd, ...
        'm2_ant', m2_ant, 'm2_ant_nd', m2_ant_nd, 'm2_ref', m2_ref, 'm2_ref_nd', m2_ref_nd, ...
        'm3_ant', m3_ant, 'm3_ant_nd', m3_ant_nd, 'm3_ref', m3_ref, 'm3_ref_nd', m3_ref_nd, ...
        'm4_ant', m4_ant, 'm4_ant_nd', m4_ant_nd, 'm4_ref', m4_ref, 'm4_ref_nd', m4_ref_nd, ...
        't3_ant', t3_ant, 't3_ant_nd', t3_ant_nd, 't3_ref', t3_ref, 't3_ref_nd', t3_ref_nd, ...
        't4_ant', t4_ant, 't4_ant_nd', t4_ant_nd, 't4_ref', t4_ref, 't4_ref_nd', t4_ref_nd);
end
if process_SB
    SubMomGroup = struct('m1_16_ant', m1_16_ant, 'm1_16_ant_nd', m1_16_ant_nd, 'm1_16_ref', m1_16_ref, 'm1_16_ref_nd', m1_16_ref_nd, ...
        'm2_16_ant', m2_16_ant, 'm2_16_ant_nd', m2_16_ant_nd, 'm2_16_ref', m2_16_ref, 'm2_16_ref_nd', m2_16_ref_nd, ...
        'm3_16_ant', m3_16_ant, 'm3_16_ant_nd', m3_16_ant_nd, 'm3_16_ref', m3_16_ref, 'm3_16_ref_nd', m3_16_ref_nd, ...
        'm4_16_ant', m4_16_ant, 'm4_16_ant_nd', m4_16_ant_nd, 'm4_16_ref', m4_16_ref, 'm4_16_ref_nd', m4_16_ref_nd, ...
        't3_16_ant', t3_16_ant, 't3_16_ant_nd', t3_16_ant_nd, 't3_16_ref', t3_16_ref, 't3_16_ref_nd', t3_16_ref_nd, ...
        't4_16_ant', t4_16_ant, 't4_16_ant_nd', t4_16_ant_nd, 't4_16_ref', t4_16_ref, 't4_16_ref_nd', t4_16_ref_nd);
end


elapsedtime = toc;

disp(['Execution Time: ' num2str(elapsedtime) ' s']);

[filepath, filename, ~] = fileparts(filename);

if process_FB && process_SB
    save(fullfile(filepath, strcat(filename, '_FB_SB', '.mat')),'FullMomGroup', 'SubMomGroup', 'time', 'pulses', 'cycles', 'mode', 'pcm', '-v7.3');
elseif process_FB
    save(fullfile(filepath, strcat(filename, '_FB', '.mat')),'FullMomGroup', 'time', 'pulses', 'cycles', 'mode', 'pcm');
else
    save(fullfile(filepath, strcat(filename, '_SB', '.mat')),'SubMomGroup', 'time', 'pulses', 'cycles', 'mode', 'pcm', '-v7.3');
end

