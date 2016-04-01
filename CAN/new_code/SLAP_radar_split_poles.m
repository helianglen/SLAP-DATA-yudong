function [time_out, radar_hh, radar_hv, radar_vh, radar_vv] = SLAP_radar_split_poles(time, radar_h, radar_v, pcm)

[radar_hh, radar_hv, radar_vh, radar_vv] = deal(nan(size(radar_h)));

% determine which data corresponds to co-pol and which corresponds to
% cross-pol using the PCM vector
% if the PCM has a 0 in its 4s bit, it's H transmitting
% if it has a 1 in its 4s bit, it's V transmitting
v_tx = bitget(pcm, 3) == 1;

radar_hh(~v_tx,:) = radar_h(~v_tx,:);
radar_hv(v_tx,:) = radar_h(v_tx,:);

radar_vh(~v_tx,:) = radar_v(~v_tx,:);
radar_vv(v_tx,:) = radar_v(v_tx,:);

time_out = time;