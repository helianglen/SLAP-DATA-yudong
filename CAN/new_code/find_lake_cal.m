%% Find and print lake calibration data
%  Find lake cal period from roll angle readings, 
%  then find Tbs from buoy data table
%    

clear

% min time-span of cal data to use -- to get rid of airplane quick turn-arounds  
cal_dt = 30;   % secs 

%sample interval: 1/15 sec
dt = 1/15.0;  

dw = cal_dt / dt;    %anything narrower than this is not lake flyover 


%lower limit: sometimes the time range went berserk 
time0 = datenum('2015-Oct-1 00:00:00', 'yyyy-mmm-dd HH:MM:SS'); 

% control plot range 
countLow = 2.5e6; 
countHigh = 4e6; 

% read in lake Tb data 
tdata = readtable('/home/ytian/smos/SLAP-DATA-yudong/CAN/Lake_buoys/buoy1_all_data.csv', ...
                  'ReadVariableNames', false); 

buoy_time = datenum(tdata{:, 1}, 'mm/dd/yyyy HH:MM');   
buoy_temp = tdata{:, 3}; 

%to plot
plot(buoy_time, buoy_temp)
datetick('x', 6)
print('lake_tb-vs-time', '-dpng');

filesRad = dir('RAD*FB_v4.mat');
files_to_process = 1:length(filesRad);

timeAll = [];
h2antAll = [];
v2antAll = [];
rollAll = [];

for i = files_to_process
   % time data from FB.mat and h2ant and v2ant from FB_m2data.mat
   load(filesRad(i).name)                       % RADTELEM_20151108T212000_FB_v4.mat
   fdate=filesRad(i).name(10:17);

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

% find lake flyover from roll angle data, and find the closest match in lake Tb data

r40=find( (rollAll > -43 & rollAll < -37) | (rollAll > 37 & rollAll < 43) );

% get a 0/1 flag vector 
roll01 = rollAll*0; 
roll01(r40) = 1; 

[pks, locs, ws, ps] = findpeaks(roll01); 

for i=1:length(ws) 
  if (ws(i) < dw ) 
   %wipe out this narrow step 
   wipe=locs(i):(locs(i)+ws(i)-1);  
   roll01(wipe) = 0; 
  end 
end 

cal=find(roll01);   % indices of usable cal data 
lake_time=timeAll(cal); 
lake_h2=h2antAll(cal); 
lake_v2=v2antAll(cal); 
lake_roll=rollAll(cal); 
lake_tb=interp1(buoy_time, buoy_temp, lake_time);

 x1= min(lake_time); 
 x2= max(lake_time); 

  figure
  subplot(3, 1, 1)
  plot(timeAll, v2antAll, 'k') 
  hold on 
  plot(lake_time, lake_v2, 'r', 'LineWidth', 5)
  axis([-Inf Inf countLow countHigh])
  title(['Time series of V-pol counts ' fdate])
  datetick('x', 15)

  subplot(3, 1, 2)
  plot(timeAll, h2antAll, 'k')
  hold on
  plot(lake_time, lake_h2, 'r', 'LineWidth', 5)
  axis([-Inf Inf countLow countHigh])
  title(['Time series of H-pol counts ' fdate])
  datetick('x', 15)

  subplot(3, 1, 3)
  plot(timeAll, rollAll, 'k') 
  hold on
  plot(lake_time, lake_roll, 'r', 'LineWidth', 5)
  axis([-Inf Inf -50 50])
  title('Time series of roll angels (deg)')
  datetick('x', 15)

  print(['lake_cal_counts_ts-' fdate], '-dpng');


  %print the cal data 

  cal_date = datestr(mean(lake_time), 'yyyy-mmm-dd HH:MM:SS');

  %disp('----------------- after outlier removal ---------------------------------------------------------------')
  fprintf('%s, %7.1f, %10.0f, %10.0f, %10.0f, %10.0f, %10.0f', ...
      cal_date, [mean(lake_tb), mean(lake_v2), mean(lake_h2), std(lake_v2), std(lake_h2), length(lake_v2)]);
  fprintf('\n');  

  %quit 
