function time = SLAP_fix_timestamps(time, pulses, cycles, sample_rate)
%SLAP_fix_timestamps Fix the timestamps finds corrupted timestamps and
%fills in the data assuming a constant sample rate.
%   [time] = SLAP_fix_timestamps(time, pulses, cycles, sample_rate) finds
%   corrupted timestamps between a missed IRIG rising edge and a pulse
%   counter reset (double 8 ms pulse) and corrects them. It also fills in
%   any NaN timestamps flagged by the IRIG decoder function
%   (SLAP_IRIG_ns_since_TAI_v2).

pri_days = sample_rate/86400;

%if the cycles goes above 1e4, assume constant sample rate until pulse
%counter resets

firstnonnan = find(~isnan(time), 1, 'first');

i = firstnonnan+1;
while i <= length(time)
    if cycles(i) > 1e4 %cycles counter overflowed which means rising edge of pulse was missed
        startpulse = pulses(i); %save the pulse value to monitor for reset
        while  i <= length(time) && pulses(i) >= startpulse %until pulse resets, increment time by sample rate
            time(i) = time(i-1) + pri_days;
            i = i+1;
        end
    else
        i = i+1;
    end
    
end

%find other timestamps flagged by the IRIG decoder and assume constant
%sample rate for them as well

for i=firstnonnan+1:length(time)
    if isnan(time(i))
        time(i) = time(i-1) + pri_days;
    end
end


end

