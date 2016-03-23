clear


counts_foambox_h = 4.1119e+06;
counts_foambox_v = 3.8788e+06;
tb_foambox = 18+273.16;

% Average of sky calibration values
skycalh = 2875917.07861282;
skycalv = 2518478.90154245;
% assumed cold sky temp in K
tb_skycal = 10;


%New DAQ columns per Albert
%TempLNAH is the IMA H Amplifier, TempLNAV is the IMA V amplifier, 
%TempAntennaH is the microstrip H 9-way and TempAntennaV is the microstrip V 9-way

load RADTELEM_20151108T163001_FB_v4

x=find(~isnan(havg) & havg>0);

plot(x, havg(x), x, vavg(x))
plot(x, azavg(x)) 

keyboard 

%% Only plot data on fore or aft half of azimuthal scans. 
 % To do this, first adjust az values to make sure they range from [0, 360]
   while min(azavg) < 0
       azavg(azavg < 0) = azavg(azavg < 0) + 360;
   end
     
   while max(azavg) > 360
       azavg(azavg > 360) = azavg(azavg > 360) - 360;
   end

plot(x, azavg(x)) 

keyboard 
        
