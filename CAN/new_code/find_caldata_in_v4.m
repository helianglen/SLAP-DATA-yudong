%% Read in the foambox temperature data, find the matching period of calibration, 
%  then use a two-round method to get rid of outliers: 
%   The first round computes mean and std from the data
%   The 2nd round gets rid of outliers (out of 3-sigma range) and recompute mean and std 

clear

% time-range  of usable cal data: if foambox temperature read at t0, use radiometer data of t0+/-cal_dt 
cal_dt = 5 / (24*60);   % 5 minutes

%lower limit: sometimes the time range went berserk 
time0 = datenum('2015-Oct-1 00:00:00', 'yyyy-mmm-dd HH:MM:SS')  

% control plot range 
foamLow = 3e6; 
foamHigh = 4e6; 

% read in foambox temperature data 
tdata = readtable('/home/ytian/smos/SLAP-DATA-yudong/CAN/edited_foam_box_data.txt'); 

cal_time = datenum(tdata{:, 1}, 'yyyy-mmm-dd HH:MM:SS')  
cal_temp = tdata{:, 2} + 273.16; 

filesRad = dir('RAD*FB_v4.mat');
files_to_process = 1:length(filesRad);

timeAll = [];
h2antAll = [];
v2antAll = [];
rollAll = [];

for i = files_to_process
    % time data from FB.mat and h2ant and v2ant from FB_m2data.mat 
    load(filesRad(i).name)                       % RADTELEM_20151005T145000_FB_m2data.mat
 
    timeAll = [timeAll; timeavg]; 
    h2antAll = [h2antAll; havg]; 
    v2antAll = [v2antAll; vavg]; 
    rollAll = [rollAll; rollavg]; 

end
  
tgood  = find(timeAll > time0); 
timeAll = timeAll(tgood);  
h2antAll = h2antAll(tgood);  
v2antAll = v2antAll(tgood);  

for j = 1:length(cal_time) 
  t_index = cal_time(j);  
  if (t_index < min(timeAll) | t_index > max(timeAll) )  % not this time 
     continue 
  end 
  % within time range
  strcat( {datestr(min(timeAll),'mm-dd HH:MM '), datestr(t_index, 'mm-dd HH:MM '), datestr(max(timeAll), 'mm-dd HH:MM ') }) 
  cal_ind = find( timeAll > (t_index - cal_dt) & timeAll < (t_index + cal_dt) ) ; 
  h2foam=h2antAll(cal_ind);  
  v2foam=v2antAll(cal_ind);  
  tfoam = timeAll(cal_ind); 

  figure
  subplot(2, 1, 1) 
  plot(timeAll, v2antAll, 'b')
  axis([min(timeAll) max(timeAll) foamLow foamHigh]) 
  title('Time series of V-pol counts')
  datetick('x', 15)
  hold on 
  plot(tfoam, v2foam, 'r') 

  subplot(2, 1, 2) 
  plot(timeAll, h2antAll, 'g')
  axis([min(timeAll) max(timeAll) foamLow foamHigh]) 
  title('Time series of H-pol counts')
  datetick('x', 15)
  hold on 
  plot(tfoam, h2foam, 'r') 

  plotname = strcat('v4_cal_ts_', datestr(t_index, 'HH_MM')); 
  print(plotname, '-dpng');  

  % get rid of empty data 
  use_ind = find( ~isnan(h2foam) ); 
  h2foam = h2foam(use_ind); 
  v2foam = v2foam(use_ind); 

  cal_date = datestr(t_index, 'yyyy-mmm-dd HH:MM:SS'); 
  disp('Date_time            ref_temp  mean(v2foam) mean(h2foam) std(v2foam) std(h2foam)     N') 
  disp('----------------- before outlier removal ---------------------------------------------------------------')
  fprintf('%s, %7.1f, %10.0f, %10.0f, %10.0f, %10.0f, %10.0f', ...
          cal_date, [cal_temp(j), mean(v2foam), mean(h2foam), std(v2foam),  std(h2foam), length(v2foam)]); 
  fprintf('\n');  

  %new ranges 
  new_h2foamLow = mean(h2foam) - 3.0 * std(h2foam);
  new_h2foamHigh = mean(h2foam) + 3.0 * std(h2foam);  

  new_v2foamLow = mean(v2foam) - 3.0 * std(v2foam);
  new_v2foamHigh = mean(v2foam) + 3.0 * std(v2foam);

  new_h2foam=h2foam(h2foam > new_h2foamLow & h2foam < new_h2foamHigh); 
  new_v2foam=v2foam(v2foam > new_v2foamLow & v2foam < new_v2foamHigh); 

  disp('----------------- after outlier removal ---------------------------------------------------------------')
  fprintf('%s, %7.1f, %10.0f, %10.0f, %10.0f, %10.0f, %10.0f', ...
          cal_date, [cal_temp(j), mean(new_v2foam), mean(new_h2foam), std(new_v2foam), std(new_h2foam), length(new_v2foam)]); 
  fprintf('\n\n');  

end % for j


