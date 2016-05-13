%% Read in the foambox temperature data, find the matching period of calibration, 
%  then use a two-round method to get rid of outliers: 
%   The first round computes mean and std from the data
%   The 2nd round gets rid of outliers (out of 3-sigma range) and recompute mean and std 

clear

% time-range  of usable cal data: if foambox temperature read at t0, use radiometer data of t0+/-cal_dt 
cal_dt = 5 / (24*60);   % 5 minutes

%lower limit: sometimes the time range went berserk
time0 = datenum('2015-Oct-1 00:00:00', 'yyyy-mmm-dd HH:MM:SS')

%plot range 
foamLow=2.2e6
foamHigh=4e6

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
rollAll = rollAll(tgood);

  x1= min(timeAll)
  x2= max(timeAll)

  figure
  subplot(3, 1, 1) 
  plot(timeAll, v2antAll, 'b')
  axis([x1 x2 foamLow foamHigh]) 
  title('Time series of V- and H-pol counts')
  datetick('x', 15)
  hold on 
  plot(timeAll, h2antAll, 'g')

  subplot(3, 1, 2) 
  plot(timeAll, h2antAll-v2antAll)
  axis([x1 x2 0e5 7e5]) 
  title('Time series of H - V counts')
  hold on 
  plot([x1 x2], [2e5 2e5], 'g')   % 2e5 is about the difference producing TbH-TbV=0. 
  datetick('x', 15)

  subplot(3, 1, 3) 
  plot(timeAll, rollAll)
  axis([min(timeAll) max(timeAll) -50 50]) 
  title('Time series of roll angels (deg)')
  datetick('x', 15)

  print('all_counts_ts', '-dpng');  

