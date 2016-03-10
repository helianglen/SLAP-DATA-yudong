function [time_out, H_I, H_Q, V_I, V_Q] = SLAP_FB_reshape(time, FB_data)

h_i_temp = squeeze(FB_data(:,1,:));
h_q_temp = squeeze(FB_data(:,2,:));
v_i_temp = squeeze(FB_data(:,3,:));
v_q_temp = squeeze(FB_data(:,4,:));

H_I = reshape(h_i_temp', size(h_i_temp, 1)*size(h_i_temp, 2), 1);
H_Q = reshape(h_q_temp', size(h_q_temp, 1)*size(h_q_temp, 2), 1);
V_I = reshape(v_i_temp', size(v_i_temp, 1)*size(v_i_temp, 2), 1);
V_Q = reshape(v_q_temp', size(v_q_temp, 1)*size(v_q_temp, 2), 1);

sample_period = 500e-6; %sec
sample_period_days = sample_period/86400;
time_temp = [time, time + sample_period_days, time + sample_period_days*2, time + sample_period_days*3];

time_out = reshape(time_temp', size(time_temp, 1) * size(time_temp, 2), 1);