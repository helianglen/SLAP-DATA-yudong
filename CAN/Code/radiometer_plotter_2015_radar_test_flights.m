% clear
% for naming of KML file
KML_file_name = ['Lines_3-4']; % land/water boundary

% determine flight date based on date in folder name
flight_directory = pwd;
flight_day = str2double((flight_directory(10:11)));
flight_month = str2double((flight_directory(8:9)));

% flag to decide whether to plot all of the scan or just half, specifically
% the fore or aft halfs
flagHalf = 1;

% flag to process data for fore or aft half-scan
fore = 1;

pol = 1; % 1 == hpol and 0 == vpol

% flag to plot h or v pol data on Google Earth
GoogleEarth  = 1;

% flag to choose whether to perform calculations to put data in tabular form
data_accum = 0;
% What to divide the footprint by to make each pixel on Google Earth.
% For visualization purposes only.
scaler = 4;
% set cal box values and time frame of interest based on the flight day
if flight_day == 12 && flight_month == 2
    Box_Cal_H_pol = 3.91e6; %+/- 0.07e6
    Box_Cal_V_pol = 3.98e6; %+/- 0.09e6
    Box_Cal_Temp = 8 + 273.16; % in Kelvin
    timeRadInit = datenum('2/12/2015 17:28:34');
    timeRadFinal = datenum('2/12/2015 17:55:55');
elseif flight_day == 25 && flight_month == 2
    Box_Cal_H_pol = 3.86e6; %+/- 0.09e6
    Box_Cal_V_pol = 4.07e6; %+/- 0.10e6
    Box_Cal_Temp = 13.4 + 273.16; % in Kelvin
    % test point 11 line 1
    % timeRadInit = datenum('2/25/2015 21:40');
    % timeRadFinal = datenum('2/25/2015 21:44');
    % % test point 11 line 2
    % timeRadInit = datenum('2/25/2015 21:46:15');
    % timeRadFinal = datenum('2/25/2015 21:48:34');
    % test point 11 line 3
    timeRadInit = datenum('2/25/2015 21:50');
    timeRadFinal = datenum('2/25/2015 21:53:31');
    % test point 11 line 4
    timeRadInit = datenum('2/25/2015 21:55:37');
    timeRadFinal = datenum('2/25/2015 21:57:57');
    % test point 12 line 1-4
    timeRadInit = datenum('2/25/2015 22:00:44');
    timeRadFinal = datenum('2/25/2015 22:15:37');
    %     % test point 14 line 1-9
    %     timeRadInit = datenum('2/25/2015 22:25:58');
    %     timeRadFinal = datenum('2/25/2015 23:06:25');
elseif flight_day == 1 && flight_month == 4
    Box_Cal_H_pol = 3.7890e+06;
    Box_Cal_V_pol = 3.4790e+06;
    Box_Cal_Temp = 15.2 + 273.16; % in Kelvin
    timeRadInit = datenum('4/1/2015 20:09:35');
    timeRadFinal = datenum('4/1/2015 20:17:02');
   
end

load radiometer_6deg_avg

% only get GPS receiver data for time period of interest
indRad = find(timeavg >= timeRadInit & timeavg < timeRadFinal );

latavg = latavg(indRad);
lonavg = lonavg(indRad);
altavg = altavg(indRad);
trkavg = trkavg(indRad);
hdgavg = hdgavg(indRad);
rollavg = rollavg(indRad);
timeavg = timeavg(indRad);
azavg = azavg(indRad);
havg = havg(indRad);
vavg = vavg(indRad);

% flag to calculate or not calculate SM
soil = 0;
% flag to plot forest data or to leave it out (for SM purposes)
noforest = 0;

% value used to remove high SM values on water
maxSM = 1;

% specify max and min values for color bar range in Google Earth KML file
if ~soil
    minBin=200;
    maxBin=300;
else
    minBin=0;
    maxBin=maxSM;
end

% specify the attributes of the current data set being worked on to append
% to filenames of output
if flagHalf
    if fore
        name = '_fore';
    else
        name = '_aft';
    end
end

if pol
    poltype = '_h';
else
    poltype = '_v';
end

% flag to decide whether to plot data during turns
flagRoll = 1;

% Aircraft roll angle limits if removing data during turns
rollmin = -10;
rollmax = 15;

% sky cal data from Jan 22, 2015 at GSFC Softball Field Complex, MD
Sky_Cal_H_pol = 2.76e6; %+/- 0.07e6
Sky_Cal_V_pol = 2.715e6; %+/- 0.065e6
Sky_Cal_Temp = 10; % in Kelvin

% initiate variables
[SM_total, az_total, emiss_surf_total, clay_total, emissinit_total, temp_h_switch, temp_v_switch, gamma_total, alt_total roll_total trk_total yaw_total hdg_total tau_total emiss_soil_total emiss_h_total  Lat_total, ndvi_total, tbh_total, tbv_total, Lon_total, ndvi_total, lcclass_total, lst_total, time_total] = deal([]);

% daq temps imported and combined into one file for entire flight in
% CombineAncillaryData.m
load daq_all

% pull off temperatures of interest
timeDaq = daq_all(:,1);
IMA_H_SWITCH = daq_all(:,9);
IMA_V_SWITCH = daq_all(:,6);
TempDiplexerH = daq_all(:, 4);
TempDiplexerV = daq_all(:, 5);
TempLNAH = daq_all(:,8);
TempLNAV = daq_all(:,7);
TempAntennaH = daq_all(:,21);
TempAntennaV = daq_all(:,22);

%% Only plot data on fore or aft half of azimuthal scans.
% To do this, first adjust az values to make sure they range from [0, 360]

azavg(azavg>360) = azavg(azavg>360) - 360;

% take into account the tracking angle so the proper first half of the
% scan can be obtained

% initialize flag_half_scan as zeros
flag_half_scan = zeros(length(trkavg),1);

for k = 1:length(trkavg)
    % get 180 degrees of scan width around track angle to get front half
    scan_lower = trkavg(k) + 90;
    scan_upper = trkavg(k) - 90;
    %     scan_lower(scan_lower>360) = scan_lower(scan_lower>360) - 360;
    %     scan_upper(scan_upper<0) = scan_upper(scan_upper<0) + 360;
    
    % get indices of az angles within the 180 degrees of the front half
    flagf1 = azavg(k)>scan_upper;
    flagf2 = azavg(k)<=scan_lower;
    
    if trkavg(k) < 180;
        
        if (flagf1 && flagf2)
            flag_half_scan(k) = k;
            % this situation is the most clear cut where the az angle falls in
            % between two values of scan_upper and scan_lower within [0 360] deg
        elseif scan_upper < 0
            % this situation requires more adjustment, where scan_upper is
            % negative so you must add 360 deg to it to make it within the
            % range of the az angles, then determine if the az angle is
            % higher than the new scan_upper, which means its in the fore HS
            if azavg(k) > (scan_upper + 360)
                flag_half_scan(k) = k;
            end
            % otherwise the az angle is in the aft half scan so don't
            % include its index in the fore half scan indices
        end
    else % trkavg > 180 deg
        % this situation is the most clear cut where the az angle falls in
        % between two values of scan_upper and scan_lower within [0 360] deg
        if (flagf1 && flagf2)
            flag_half_scan(k) = k;
            % this situation requires more adjustment, where scan_lower is >
            % 360 so it must be subtracted by 360 to fit within the
            % range of the az angles, then determine if the az angle is
            % lower than the new scan_lower, which means its in the fore HS.
        elseif scan_lower > 360
            if azavg(k) < (scan_lower - 360)
                flag_half_scan(k) = k;
            end
            % otherwise the az angle is in the aft half scan so don't
            % include its index in the fore half scan indices
        end
    end
end

% turn flag_half_scan into a vector of logicals where ones are indices of
% fore half scan data points
flag_half_scan = flag_half_scan>0;


if flagHalf
    if ~fore
        % get indices of aft half scan as the opposite of the indices of the fore
        % half scan if the flag for the fore half scan is turned off
        flag_half_scan = ~flag_half_scan;
    end
else
    % to plot fore and aft half scans together, get all logical ones
    flag_half_scan = find(latavg>0);
end

% get geolocation data for fore or aft half scans using flag_half_scan
trk_half_scan = trkavg(flag_half_scan);
lon_half_scan = lonavg(flag_half_scan);
lat_half_scan = latavg(flag_half_scan);
roll_half_scan = rollavg(flag_half_scan);
v_half_scan = vavg(flag_half_scan);
h_half_scan = havg(flag_half_scan);
alt_half_scan = altavg(flag_half_scan);
time_half_scan = timeavg(flag_half_scan);
az_half_scan = azavg(flag_half_scan);
hdg_half_scan = hdgavg(flag_half_scan);


% use cold sky and foam box observations to do a two-point
% calibration between raw counts values and Tb in K
tbcalh = interp1([Sky_Cal_H_pol Box_Cal_H_pol], [Sky_Cal_Temp,  Box_Cal_Temp], h_half_scan, 'linear', 'extrap');
tbcalv = interp1([Sky_Cal_V_pol Box_Cal_V_pol], [Sky_Cal_Temp,  Box_Cal_Temp], v_half_scan, 'linear', 'extrap');

% apply correction for galaxy and cosmic microwave background radiation
tbcalh = tbcalh - 1;
tbcalv = tbcalv - 1;

% apply net sidelobe correction
tbcalh = tbcalh*(1 - 0.115) - 0.1;
tbcalv = tbcalv*(1 - 0.082) - 0.1;


% administer roll filter
if flagRoll
    indRoll = find(roll_half_scan < rollmax & roll_half_scan > rollmin);
else
    % include all data when turning using this flag
    indRoll = logical(ones(length(h_half_scan),1));
end

roll_half_scan = roll_half_scan(indRoll);

tbh_half_scan = tbcalh(indRoll);
tbv_half_scan = tbcalv(indRoll);
time_half_scan= time_half_scan(indRoll);
trk_half_scan = trk_half_scan(indRoll);
lon_half_scan = lon_half_scan(indRoll);
lat_half_scan = lat_half_scan(indRoll);
alt_half_scan = alt_half_scan(indRoll);
az_half_scan = az_half_scan(indRoll);
hdg_half_scan = hdg_half_scan(indRoll);

% remove heading angle from az scan angle before az scan angle is
% written into the data file
az_wo_hdg = az_half_scan - hdg_half_scan;
az_wo_hdg(az_wo_hdg < 0) = az_wo_hdg(az_wo_hdg < 0) + 360;

% loss coefficients from Ed's documentation of Tb corrections
loss1 = 0.3543;
loss2 = 0.0450;
loss3 = 0.0689;
loss4 = 0.045;

% apply Tb corrections for H-Pol
Tb = tbh_half_scan;
tb_otherpol = tbv_half_scan;
Tb = Tb*(1 - 0.01) + tb_otherpol*0.01;

temp_dip_interp    = interp1(timeDaq,TempDiplexerH, time_half_scan);
temp_lna_interp    = interp1(timeDaq,TempLNAH, time_half_scan);
temp_ant_interp    = interp1(timeDaq,TempAntennaH, time_half_scan);
temp_switch_interp = interp1(timeDaq,IMA_H_SWITCH, time_half_scan);

temp_dip_interp(isnan(temp_dip_interp)) = nanmean(temp_dip_interp);
temp_lna_interp(isnan(temp_lna_interp)) = nanmean(temp_lna_interp);
temp_ant_interp(isnan(temp_ant_interp)) = nanmean(temp_ant_interp);
temp_switch_interp(isnan(temp_switch_interp)) = nanmean(temp_switch_interp);
temp_switch_interp_h = temp_switch_interp;

Tb = Tb*(1- loss1) + (loss1)*( temp_ant_interp+273.16);
Tb = Tb*(1- loss2) + (loss2)*( temp_dip_interp+273.16);
Tb = Tb*(1- loss3) + (loss3)*( temp_dip_interp+273.16);
Tbh = Tb*(1- loss4) + (loss4)*( temp_lna_interp+273.16);

% apply Tb corrections for V-Pol
Tb = tbv_half_scan;
tb_otherpol = tbh_half_scan;
Tb = Tb*(1 - 0.01) + tb_otherpol*0.01;

temp_dip_interp    = interp1(timeDaq,TempDiplexerV, time_half_scan);
temp_lna_interp    = interp1(timeDaq,TempLNAV, time_half_scan);
temp_ant_interp    = interp1(timeDaq,TempAntennaV, time_half_scan);
temp_switch_interp = interp1(timeDaq,IMA_V_SWITCH, time_half_scan);

temp_dip_interp(isnan(temp_dip_interp)) = nanmean(temp_dip_interp);
temp_lna_interp(isnan(temp_lna_interp)) = nanmean(temp_lna_interp);
temp_ant_interp(isnan(temp_ant_interp)) = nanmean(temp_ant_interp);
temp_switch_interp(isnan(temp_switch_interp)) = nanmean(temp_switch_interp);
temp_switch_interp_v = temp_switch_interp;

Tb = Tb*(1- loss1) + (loss1)*( temp_ant_interp+273.16);
Tb = Tb*(1- loss2) + (loss2)*( temp_dip_interp+273.16);
Tb = Tb*(1- loss3) + (loss3)*( temp_dip_interp+273.16);
Tbv = Tb*(1- loss4) + (loss4)*( temp_lna_interp+273.16);


% accumulate all data from 10 minute data sets into one vector for each variable
tbh_total = [tbh_total; Tbh];
tbv_total = [tbv_total; Tbv];
time_total = [time_total; datevec(time_half_scan)];
Lat_total = [Lat_total; lat_half_scan];
Lon_total = [Lon_total; lon_half_scan];
alt_total = [alt_total; alt_half_scan];
roll_total = [roll_total; roll_half_scan ];
trk_total = [trk_total; trk_half_scan ];
hdg_total = [hdg_total; hdg_half_scan ];
az_total = [az_total; az_wo_hdg];
yaw_total = [yaw_total; hdg_total - trk_total];

temp_h_switch = [temp_h_switch; temp_switch_interp_h];
temp_v_switch = [temp_v_switch; temp_switch_interp_v];

if soil
    [lcclass, lcvalue, ndvival, lst, lst_time, clay] = findSMAncillaryValues(lon_half_scan,lat_half_scan, flight_day);
    lcclass_total = [lc_class_total; lcclass];
    ndvi_total = [ndvi_total; ndvival];
    lst_total = [lst_total; lst];
    clay_total = [clay_total; clay];
end

% only accumlate tabular data if the flag to do so, data_accum, is set
if  data_accum
    
    twoSigmah = 2* nanstd(temp_h_switch);
    twoSigmav = 2* nanstd(temp_v_switch);

    switchflagh = abs(temp_h_switch - mean(temp_h_switch)) > twoSigmah;
    switchflagv = abs(temp_v_switch - mean(temp_v_switch)) > twoSigmav;

    % flag for abnormally low values (below 3 sigma from the mean of
    % the h and v pol data)
    
    flagLowValH = tbh_total < (nanmean(tbh_total) - 3 * nanstd(tbh_total));
    flagLowValV = tbv_total < (nanmean(tbv_total) - 3 * nanstd(tbv_total));
    
    % flag for RFI above 3 sigma from the mean of the h and v pol data
    
    tbhnorm = tbh_total(~flagLowValH);
    tbh_total_rfi_threshold = nanmean(tbhnorm) + 3* nanstd(tbhnorm);
    flagRFI_hpol = tbh_total>tbh_total_rfi_threshold;
    
    tbvnorm = tbv_total(~flagLowValV);
    tbv_total_rfi_threshold = nanmean(tbvnorm) + 3* nanstd(tbvnorm);
    flagRFI_vpol = tbv_total>tbv_total_rfi_threshold;
    
    % set flag to determine if SLAP is spinning or not based on the
    % difference in consecutive azimuth scan angles being less than 4 deg
    scan_motor_diff = abs(diff(az_total));
    scan_motor_spinning = scan_motor_diff > 4;
    scan_motor_spinning = [scan_motor_spinning; logical(0)];
    
    % determine if data point is in the fore or aft half scan 
    % (N/A to stare data)
    flagHS = az_total > 270 | az_total < 90;
    
    % combine all data into one matrix
    alldata = [ time_total, Lat_total, Lon_total, alt_total, tbh_total, tbv_total,  az_total, roll_total, trk_total, hdg_total, switchflagh, switchflagv, flagLowValH, flagLowValV, flagRFI_hpol, flagRFI_vpol, scan_motor_spinning, flagHS];
    
    % remove data lines with all NaN values
    indnan = isnan(Lat_total);
    alldata = alldata(~indnan,:);
    
    % save as csv file
    csvwrite([num2str(flightdate(4:11)) '  Observations.csv'] , alldata)
end

% choose polarization of Tb data for Soil Moisture processing and/or Google Earth
% image generation
if pol
    Tb = tbh_total;
else
    Tb = tbv_total;
end

if soil
    [SM , emiss_init, emiss_surf , emiss_soil , emiss_h, tau, gamma, lstnew ] = deal(zeros(length(Tb),1));
    lstnonzero = lst_total(lst_total>0);
    lst_total(lst_total==0) = mean(lstnonzero);
    
    for k = 1:length(Tb)
        [SM(k),emiss_init(k), emiss_surf(k), emiss_soil(k), emiss_h(k), tau(k), gamma(k)] = Tb2SM(Tb(k), lst_total(k), clay_total(k)*100, 0:0.005:1, ndvi_total(k), lc_class_total(k));
    end
    
    if noforest
        SM(lc_class_total<6) = nan;
    end
    Tb = SM;
end
% plot data as an image on Google Earth
if GoogleEarth
    
    lat_filtered = Lat_total;
    lon_filtered = Lon_total;
    alt_filtered = alt_total;
    roll_filtered = roll_total;
    % SLAP's elevation angle is a function of the plane's roll angle. When
    % the plane is flying straight and level, the elevation, or incidence,
    % angle is 40 deg.
    elev = 40;
    AGL = alt_filtered;
    
    % determine footprint size based on altitude
    slant = AGL./ cosd(elev);
    
    if pol
        beam_width = 18.8;
    else
        beam_width = 18.1;
    end
    footprintm = (beam_width./180.*3.14.*slant)./cosd(elev);
    r_earth = 6378e3;

    % convert footprint radius in meters to degrees longitude
    footprint_major = footprintm / (pi*r_earth) * 180/scaler;
    footprint_minor = footprint_major;
    
    
    % generate KML file to plot on Google Earth with lat/lon and Tb/SM data
    
    fida=fopen([KML_file_name, name, poltype,  '.kml'],'w');
    
    line1=['<?xml version="1.0" encoding="UTF-8"?>'];
    line2=['<kml xmlns="http://earth.google.com/kml/2.0">'];
    line3=['<Document>'];
    line4=[strcat('<name>',[KML_file_name, name, poltype, '.kml'],'</name>')];
    line6=['<Style id="air"><icon><href>root://icons/palette-4.png?x=160&amp;y=0&amp;w=32&amp;h=32</href></icon>'];
    line7=['<LineStyle><color>FF0000FF</color><width>2.0</width></LineStyle></Style>'];
    line8=['<Style id="gnd"><icon><href>root://icons/palette-4.png?x=160&amp;y=0&amp;w=32&amp;h=32</href></icon>'];
    line9=['<LineStyle><color>FFBBFF00</color><width>2.0</width></LineStyle></Style>'];
    line10=['<IconStyle><color>FFFFFFFF</color><scale>0.7</scale></IconStyle>'];
    line11=['<LabelStyle><color>FFFFFFFF</color><scale>0.7</scale></LabelStyle></Style>'];
    line_folder=['</Folder>'];
    line_last=['</Document></kml>'];
    
    
    cls=['FF0000AA';'FF0000D4';'FF0000FF';'FF002AFF';'FF0055FF';'FF0080FF';'FF00AAFF'; ...
        'FF00D4FF';'FF00FFFF';'FF2AFFD4';'FF54FFAA';'FF80FF80';'FFAAFF55';'FFD4FF2A'; ...
        'FFFFFF00';'FFFFD400';'FFFFAA00';'FFFF8000';'FFFF5500';'FFFF2A00';'FFFF0000';'50000000';'50000000'];
    
    if ~soil
        cls =     flipud(cls);
    end
    
    fprintf(fida,'%s\r%s\r%s\r%s\r%s\r%s\r%s\r%s\r', ...
        line1,line2,line3,line4,line6,line7,line8,line9);
    
    %-----------------------------------
    % bin into value ranges (20 colours and 2 outside of range)
    %-----------------------------------
    
    deltaB=(maxBin-minBin)./20;
    f0=find(Tb<minBin); % | isnan(Tb));
    f1=find(Tb>=minBin & Tb<minBin+deltaB*1);
    f2=find(Tb>=minBin+deltaB*1 & Tb<minBin+deltaB*2);
    f3=find(Tb>=minBin+deltaB*2 & Tb<minBin+deltaB*3);
    f4=find(Tb>=minBin+deltaB*3 & Tb<minBin+deltaB*4);
    f5=find(Tb>=minBin+deltaB*4 & Tb<minBin+deltaB*5);
    f6=find(Tb>=minBin+deltaB*5 & Tb<minBin+deltaB*6);
    f7=find(Tb>=minBin+deltaB*6 & Tb<minBin+deltaB*7);
    f8=find(Tb>=minBin+deltaB*7 & Tb<minBin+deltaB*8);
    f9=find(Tb>=minBin+deltaB*8 & Tb<=minBin+deltaB*9);
    f10=find(Tb>=minBin+deltaB*9 & Tb<=minBin+deltaB*10);
    f11=find(Tb>=minBin+deltaB*10 & Tb<=minBin+deltaB*11);
    f12=find(Tb>=minBin+deltaB*11 & Tb<=minBin+deltaB*12);
    f13=find(Tb>=minBin+deltaB*12 & Tb<=minBin+deltaB*13);
    f14=find(Tb>=minBin+deltaB*13 & Tb<=minBin+deltaB*14);
    f15=find(Tb>=minBin+deltaB*14 & Tb<=minBin+deltaB*15);
    f16=find(Tb>=minBin+deltaB*15 & Tb<=minBin+deltaB*16);
    f17=find(Tb>=minBin+deltaB*16 & Tb<=minBin+deltaB*17);
    f18=find(Tb>=minBin+deltaB*17 & Tb<=minBin+deltaB*18);
    f19=find(Tb>=minBin+deltaB*18 & Tb<=minBin+deltaB*19);
    f20=find(Tb>=minBin+deltaB*19 & Tb<=minBin+deltaB*20);
    f21=find(Tb>maxBin);
    
    for j=1:22
        fprintf(fida,'%s %s %s %s %s\r','   <Placemark><name>TB Range: ',num2str((minBin+deltaB*(j-2))) ...
            ,'-',num2str((minBin+deltaB*(j-1))),'(K)</name>');
        fprintf(fida,'%s\r','       <Style>');
        fprintf(fida,'%s\r','           <LineStyle>');
        fprintf(fida,'%s %s %s\r','               <color>',cls(j,:),'</color>');
        fprintf(fida,'%s\r','               <width> 0 </width>');
        fprintf(fida,'%s\r','           </LineStyle>');
        fprintf(fida,'%s\r','           <PolyStyle>');
        fprintf(fida,'%s %s %s\r','               <color>',cls(j,:),'</color>');
        fprintf(fida,'%s\r','           </PolyStyle>');
        fprintf(fida,'%s\r','       </Style>');
        fprintf(fida,'%s\r','     <MultiGeometry>');
        fbin = eval(genvarname(['f',num2str(j-1)]));
        
        if ~isempty(fbin)
            for n=1:length(fbin)
                
                fprintf(fida,'%s\r','       <Polygon>');
                fprintf(fida,'%s\r','           <altitudeMode>relativetoground</altitudeMode>');
                fprintf(fida,'%s\r','           <outerBoundaryIs>');
                fprintf(fida,'%s\r','               <LinearRing>');
                fprintf(fida,'%s\r','                   <coordinates>');
                
                % create ellipse of observations
                for theta = 0:15:360
                    fprintf(fida,'%s%6f,%6f,%0f\r','                   ',lon_filtered(fbin(n))+footprint_minor(fbin(n))*cosd(theta),lat_filtered(fbin(n))+footprint_major(fbin(n))*sind(theta), 6000);
                end
                
                fprintf(fida,'%s\r','                   </coordinates>');
                fprintf(fida,'%s\r','               </LinearRing>');
                fprintf(fida,'%s\r','           </outerBoundaryIs>');
                fprintf(fida,'%s\r','       </Polygon>');
            end
        end
        
        fprintf(fida,'%s\r','     </MultiGeometry>');
        fprintf(fida,'%s\r','   </Placemark>');
        
    end
    
    fprintf(fida,'%s\r',line_last);
    fclose(fida);
    
end