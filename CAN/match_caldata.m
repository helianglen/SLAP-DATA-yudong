
function [lake_tb, lake_v, lake_h, foam_tb, foam_v, foam_h] = match_caldata(input_time)    

% Given a time as input, find the closest foambox and lake calibration data 

%testing 
%input_time = datenum('2015-11-08 10:00:00', 'yyyy-mmm-dd HH:MM:SS');

%foam box 
cal=readtable('/home/ytian/smos/SLAP-DATA-yudong/CAN/Cal_table.csv'); 
foam_time = datenum(cal{:, 1}, 'yyyy-mmm-dd HH:MM:SS');
tb = cal{:, 2}; 
yv = cal{:, 3}; 
yh = cal{:, 4}; 

%lake 
lcal=readtable('/home/ytian/smos/SLAP-DATA-yudong/CAN/Lake_cal_table.csv'); 
lake_time = datenum(lcal{:, 1}, 'yyyy-mmm-dd HH:MM:SS');
ltb = lcal{:, 2}; 
lyv = lcal{:, 3}; 
lyh = lcal{:, 4}; 

[tc index]=min(abs(foam_time - input_time)); 
foam_tb = tb(index); 
foam_v = yv(index); 
foam_h = yh(index); 

[tc index]=min(abs(lake_time - input_time)); 
lake_tb = ltb(index); 
lake_v = lyv(index); 
lake_h = lyh(index); 


