function [time, pulses, cycles] = SLAP_IRIG_ns_since_TAI_v2(irig_bytes, swap_endianness)

irig_bits = reshape(dec2bin(irig_bytes, 8)', 1, 128);

%check to see if the IRIG packet framing is all zeroes like it should be.
%if not, the irig frame is probably bad
if all([irig_bits(1:2) irig_bits(11:10:91)] == '0');
    
    cycles = bin2dec(irig_bits(109:128));
    
    pulses = bin2dec(irig_bits(101:108));
    
    %FPGA stores sequential bits received MSB-first, and dec2bin is also
    %big-endian, so index 1 is the first bit received in the IRIG frame
    irig_code = irig_bits(1:100);
    
    %pull out the data from the IRIG packet framing
    ns_bits = [irig_code(3:10) irig_code(12:20) irig_code(22:30) irig_code(32:40) irig_code(42:50) irig_code(52:60) irig_code(62:70) irig_code(72:73)];
    
    %bits were sent LSB-first, so flip to use the big-endian bin2dec. this also
    %makes the entire bytestream big-endian when casting the byte array
    ns_bits = reshape(flip(ns_bits), 8, 8)';
    ns_bytes = bin2dec(ns_bits);
    ns_since_TAI = double(SLAP_UINT64(ns_bytes, swap_endianness)); %massive loss of precision here. fix it?
    
    sec_since_TAI = ns_since_TAI/1e9;
    
    
    %second rollover happens some time in the 99th pulse
    %the cutoff is determined by looking at the empirical data
    %it's a number greater than 2039 and less than 2041, at last check
    if(pulses == 99 && cycles >= 2039)
        %fracsec = 0.98 + cycles/1e6;
        fracsec = (pulses-1)/100 + cycles/1e6;
    else
        fracsec = 1 + (pulses-1)/100 + cycles/1e6;
    end
    
    sec_since_TAI = sec_since_TAI + fracsec;
    days_since_TAI = sec_since_TAI/86400;
    
    TAI_epoch = 719529;
    
    time = TAI_epoch + days_since_TAI;
    
else
    time = nan;
    pulses = nan;
    cycles = nan;
end


%if the cycle count overflows, the pulse rising edge wasn't detected and
%the timestamp is invalid
if cycles > 1e4
    time = nan;
end

%invalidate timestamps in the deadband while the seconds data increments
if pulses == 99 && cycles > 2035 && cycles < 2045
    time = nan;
end




