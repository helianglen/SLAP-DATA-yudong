clear
% flag to choose whether to perform calculations to put data in tabular form
data_accum = 1;
% flag to plot h or v pol data on Google Earth
GoogleEarth  = 0;
pol = 1; % 1 == hpol and 0 == vpol
% flag to process data for fore or aft half-scan
fore = 1;
% flag to decide whether to plot all of the scan or just half, specifically
% the fore or aft halfs
flagHalf = 0;
% load RAD files 
filesRad2 = dir('RAD*FB_v4.mat');

% choose whether to plot only NE or SW facing part of flight lines for May
% 2014 iPHEX Campaign AM flights over Western NC. If the NEorSW variable is enabled,
% the plotting algorithm chooses one side or another to plot based on flag (NE = 1, SW = 0).
% NEorSW = 1;

% flag to plot forest data or to leave it out (for SM purposes)
noforest = 0;

% flag to remove streaks in h-pol Tb data
removestreaks = 1;

% flag to calculate or not calculate SM
soil = 0;

flightdatestr = pwd;
flightdate = str2double((flightdatestr(end-8:end-7)));
% specify what data set is current being worked on to append to filenames
% of output
if fore
    name = '_Fore';
else 
    name = '_Aft';
end

% which type of data used in processing, currently M2
power = '2';



% Aircraft roll angle limits if removing data during a turn
flagRoll = 0;
rollmin = -4;
rollmax = 2;

% value used to remove high SM values on water
maxSM = 1;

if exist('NEorSW')
    if NEorSW
        if pol
            poltype = '_h_NE';
        else
            poltype = '_v_NE';
        end
    else
        if pol
            poltype = '_h_SW';
        else
            poltype = '_v_SW';
        end
    end
else
    if pol
        poltype = '_h';
    else
        poltype = '_v';
    end
end

% specify max and min values for color bar range
if ~soil
    minBin=200;
    maxBin=300;
else 
    minBin=0;
    maxBin=maxSM;
end

if flightdate == 2 % May 2, 2014 flight (First science flight)
    
    % incorporate average of in-situ Jordan Lake temp measurements
    % (subtracting offset of temp measurement gauge which on this day was 0.3 C)
    temp_lake = mean([ 19.3 19 19 19 18.8 18.9 19.2 19.4 19.2]) - 0.3;
    
    % SLAP counts when observing the lake, from data when flying over the
    % lake. Corrected for all issues noted in SLAP L1B TB documentation.
    % Minimum counts value when flying over the lake area was chosen as it
    % reasonably is expected to coincide with the observation of least land contamination
    % (which would have the effect of increaing counts values)
    counts_lake_h = 3.3328e+06;
    counts_lake_v = 2.9705e+06;
    
    % SLAP average counts when observing the foam box pre-flight
    counts_foambox_h = 4.0755e+06;
    counts_foambox_v = 3.8327e+06;
    
    % foam box temp assumed as ambient temp in Langley airport hangar during calibration
    tb_foambox = 19+273.16;
    
elseif flightdate == 5 % May 5 flight (Second science flight)

    temp_lake = mean([21.7 21.2 21.2 21.4 21.3 21.1]) - 0.9;
    
    counts_lake_h = 3.2932e+06;
    counts_lake_v = 2.8763e+06;
    
    counts_foambox_h = 4.1393e+06;
    counts_foambox_v = 3.9365e+06;
    
    tb_foambox = 15+273.16;
   
elseif flightdate == 9 % May 9 flight (Third science flight)

    temp_lake = mean([23.2 22.4 21.9 22.4 21.8 22.2 21.9 22.5 22.1]);
    
    counts_lake_h = 3.4331e+06;
    counts_lake_v = 2.9311e+06;
    
    counts_foambox_h = 4.1209e+06;
    counts_foambox_v = 3.9129e+06;
    
    tb_foambox = 24+273.16;
    
elseif flightdate == 20 % May 20 flight (pre-science radar checkout flight)
        % (no lake observations on non-science flights)
    counts_foambox_h = 3.6269e+06;
    counts_foambox_v = 3.7168e+06;
    tb_foambox = 18+273.16;
    
elseif flightdate == 21 % May 21 flight (last science flight including radar)

    temp_lake = 22.00588;
    % first flight (Langley to RDU)
%     counts_lake_h = 2933558.52019094;
%     counts_lake_v = 2883552.71435995;
   
    % second_flight (RDU to Langley)
    counts_lake_h = 2.9953e+06;
    counts_lake_v = 2.8750e+06;
    
    counts_foambox_h = 4.1209e+06;
    counts_foambox_v = 3.9129e+06;
    
%     counts_foambox_h = 3728429.78034618;
%     counts_foambox_v = 3750786.90278617;
    
    tb_foambox = 20+273.16;
    
elseif flightdate == 25 % April 25th flight (pre-science checkout and water cal flight)

    counts_foambox_h = 4.1119e+06;
    counts_foambox_v = 3.8788e+06;
    tb_foambox = 18+273.16;
    
    counts_lake_h = 3.1119e+06;
    counts_lake_v = 2.8158e+06;
end

% Take in indep lake temp measurements in Celsius and convert to Tb in K.
% Inputs are temp of the lake then 1 for lake and 0 for ocean for salinity 

% lake_or_ocean = 1;
% [temp_lake_tbh, temp_lake_tbv] = ts2tb(temp_lake, lake_or_ocean);

% Average of sky calibration values from April 21st observations
skycalh = 2875917.07861282;
skycalv = 2518478.90154245;
% assumed cold sky temp in K
tb_skycal = 10;

% initiate variables
[SM_total, az_total, emiss_surf_total, clay_total, emissinit_total, temphswitch, tempvswitch, gamma_total, alt_total roll_total trk_total hdg_total tau_total emiss_soil_total emiss_h_total  Lat_total, ndvi_total, tbh_total, tbv_total, Lon_total, ndvi_total, lcclass_total, lst_total, time_total] = deal([]);
% save SMvals SM_total emissinit_total Lat_total gamma_total tau_total emiss_surf_total emiss_soil_total emiss_h_total Lon_total lcclass_total ndvi_total  lst_total  time_total

flightmonth = str2double(flightdatestr(9));

load daq_all

timeDaq = daq_all(:,1);
IMA_H_SWITCH = daq_all(:,13);
IMA_V_SWITCH = daq_all(:,10);
TempDiplexerH = daq_all(:, 8);
TempDiplexerV = daq_all(:, 9);
TempLNAH = daq_all(:,12);
TempLNAV = daq_all(:,11);
TempAntennaH = daq_all(:,25);
TempAntennaV = daq_all(:,26);

for i = 1:length(filesRad2)
    i
    clear noscan
        load(filesRad2(i).name)
%         azavg = azavg + hdgavg;
        timestr = filesRad2(i).name(19:22);
        
        %% Only plot data on fore or aft half of azimuthal scans. 
        % To do this, first adjust az values to make sure they range from [0, 360]
        while min(azavg) < 0
            azavg(azavg < 0) = azavg(azavg < 0) + 360;
        end
        
        while max(azavg) > 360
            azavg(azavg > 360) = azavg(azavg > 360) - 360;
        end
        
        % take into account the tracking angle so the proper first half of the
        % scan can be plotted on Google Earth

        % initialize flag_half_scan as zeros
        flag_half_scan = zeros(length(trkavg),1);
        
        % only separate data into half-scans when the 'noscan' flag is
        % turned off
        if exist('noscan') ~= 1 && flagHalf
        for k = 1:length(trkavg)
            % get 180 degrees of scan width around track angle to get front half
            scanright = trkavg(k) + 90;
            scanleft = trkavg(k) - 90;
        %     scanright(scanright>360) = scanright(scanright>360) - 360;
        %     scanleft(scanleft<0) = scanleft(scanleft<0) + 360;
        
            % get indices of az angles within the 180 degrees of the front half
            flagf1 = azavg(k)>scanleft;
            flagf2 = azavg(k)<=scanright;
        
            if trkavg(k) < 180;
        
                if (flagf1 && flagf2)
                    flag_half_scan(k) = k;
                % this situation is the most clear cut where the az angle falls in
                % between two values of scanleft and scanright within [0 360] deg
                elseif scanleft < 0
                    % this situation requires more adjustment, where scanleft is
                    % negative so you must add 360 deg to it to make it within the
                    % range of the az angles, then determine if the az angle is
                    % higher than the new scanleft, which means its in the fore HS
                    if azavg(k) > (scanleft + 360)
                        flag_half_scan(k) = k;
                    end
                    % otherwise the az angle is in the aft half scan so don't
                    % include its index in the fore half scan indices
                end
            else % trkavg > 180 deg
                % this situation is the most clear cut where the az angle falls in
                % between two values of scanleft and scanright within [0 360] deg
                if (flagf1 && flagf2)
                    flag_half_scan(k) = k;
                    % this situation requires more adjustment, where scanright is >
                    % 360 so it must be subtracted by 360 to fit within the
                    % range of the az angles, then determine if the az angle is
                    % lower than the new scanright, which means its in the fore HS.
                elseif scanright > 360
                    if azavg(k) < (scanright - 360)
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
        end
        
        if exist('noscan') ~= 1 && flagHalf
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
        tbcalh = interp1([skycalh counts_foambox_h], [tb_skycal,  tb_foambox], h_half_scan, 'linear', 'extrap');
        tbcalv = interp1([skycalv counts_foambox_v], [tb_skycal,  tb_foambox], v_half_scan, 'linear', 'extrap');
        
        % apply correction for galaxy and cosmic microwave background radiation
        tbcalh = tbcalh - 1;
        tbcalv = tbcalv - 1;
        
        % apply net sidelobe correction
        tbcalh = tbcalh*(1 - 0.115) - 0.1;
        tbcalv = tbcalv*(1 - 0.082) - 0.1;
        %% administer roll filter
        if flagRoll
            indRoll = find(roll_half_scan < rollmax & roll_half_scan > rollmin);
        else
            % include all data when turning using this flag
            indRoll = 1:length(h_half_scan);
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
        
        %% separate azimuthal scan angle from heading angle
        az_wo_hdg = az_half_scan - hdg_half_scan;
        az_wo_hdg(az_wo_hdg < 0) = az_wo_hdg(az_wo_hdg < 0) + 360;
        
        %% attempt to remove streaks of high h-pol counts values at specific scan angles
        if pol && removestreaks
            [tbh_half_scan, peakval] = removetbpeaks(tbh_half_scan, az_wo_hdg);
        end
        %% apply Tb corrections as specified by Ed Kim
        loss1 = 0.3543;
        loss2 = 0.0450;
        loss3 = 0.0689;
        loss4 = 0.045;
        
        % apply corrections for H-Pol
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
        temp_switch_interph = temp_switch_interp;
        
        Tb = Tb*(1- loss1) + (loss1)*( temp_ant_interp+273.16);
        Tb = Tb*(1- loss2) + (loss2)*( temp_dip_interp+273.16);
        Tb = Tb*(1- loss3) + (loss3)*( temp_dip_interp+273.16);
        Tbh = Tb*(1- loss4) + (loss4)*( temp_lna_interp+273.16);
        
        % do corrections for V-Pol
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
        temp_switch_interpv = temp_switch_interp;
        
        Tb = Tb*(1- loss1) + (loss1)*( temp_ant_interp+273.16);
        Tb = Tb*(1- loss2) + (loss2)*( temp_dip_interp+273.16);
        Tb = Tb*(1- loss3) + (loss3)*( temp_dip_interp+273.16);
        Tbv = Tb*(1- loss4) + (loss4)*( temp_lna_interp+273.16);
        
        
        % accumulate all 10 minute data sets into one vector for each variable
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
        
        temphswitch = [temphswitch; temp_switch_interph];
        tempvswitch = [tempvswitch; temp_switch_interpv];
        
     
        % if converting Tb to soil moisture, these ancillary properties
        % must be found for each observation location
        if soil
            [lcclass, lcvalue, ndvival, lst, lst_time, clay] = findSMAncillaryValues(lon_half_scan,lat_half_scan, flightdate);
            lcclass_total = [lcclass_total; lcclass];
            ndvi_total = [ndvi_total; ndvival];
            lst_total = [lst_total; lst];
            clay_total = [clay_total; clay];
        end
        
end

% only write tabular data if the flag is set
if  data_accum
    
    twoSigmah = 2 * nanstd(temphswitch);
    twoSigmav = 2 * nanstd(tempvswitch);
    
    % a flag is triggered when the H and V Pol Switch Temps are outside
    % two standard deviations from the mean
    switchflagh = abs(temphswitch - mean(temphswitch)) > twoSigmah;
    switchflagv = abs(tempvswitch - mean(tempvswitch)) > twoSigmav;

    % flag for abnormally low values is triggered when the H and V Pol Tb
    % values are below three standard deviations from the mean
    
    flagLowValH = tbh_total< (nanmean(tbh_total) - 3 * nanstd(tbh_total));
    flagLowValV = tbv_total< (nanmean(tbv_total) - 3 * nanstd(tbv_total));
    
    % flag for RFI is triggered when the H and V Pol Tb
    % values are above three standard deviations from the mean
    
    tbh_norm = tbh_total(~flagLowValH);
    tbh_total_rfi_threshold = nanmean(tbh_norm) + 3* nanstd(tbh_norm);
    flagRFI_hpol = tbh_total>tbh_total_rfi_threshold;
    
    tbv_norm = tbv_total(~flagLowValV);
    tbv_total_rfi_threshold = nanmean(tbv_norm) + 3* nanstd(tbv_norm);
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
    csvwrite(['May ', flightdatestr(10:11), ' AM Observations.csv'] , alldata)
end

% determine which polarization to use for conversion to soil moisture and
% plotting with Google Earth based on flag
if pol
    Tb = tbh_total;
else
    Tb = tbv_total;
end

if soil
    % determine soil moisture values using the single-channel algorithm in the
    % SMAP L2_P Radiometer ATBD Document.
    [SM , emiss_init, emiss_surf , emiss_soil , emiss_h, tau, gamma, lstnew ] = deal(zeros(length(Tb),1));
    lstnonzero = lst_total(lst_total>0);
    lst_total(lst_total==0) = mean(lstnonzero);

    for k = 1:length(Tb)
        [SM(k),emiss_init(k), emiss_surf(k), emiss_soil(k), emiss_h(k), tau(k), gamma(k)] = Tb2SM(Tb(k), lst_total(k), clay_total(k)*100, 0:0.005:1, ndvi_total(k), lcclass_total(k));
    end
    
    if noforest
        SM(lcclass_total<6) = nan;
    end
    Tb = SM;    
end

if GoogleEarth
    % if the variable NEorSW exists, only plot NE or SW facing flight lines
    % based on flag (NE = 1, SW = 0)
    if exist('NEorSW')
        
        if NEorSW
            flagdir = find(trk_total > 50 & trk_total < 60);         
        else
            flagdir = find(trk_total > 230 & trk_total < 240);
        end
    else % create a flag of all ones to plot all data
        flagdir = Lat_total > 0;
    end
    
    Tb = Tb(flagdir);
    lat_filtered = Lat_total(flagdir);
    lon_filtered = Lon_total(flagdir);
    alt_filtered = alt_total(flagdir);
    roll_filtered = roll_total(flagdir);
    % elevation angle is a function of the plane's roll ange
    elev = 40 + roll_filtered;
    % AGL = Altitude above ground level
    AGL = alt_filtered;
    % determine footprint size based on altitude
    slant = AGL./ cosd(elev);
    % beam width for radiometer specified in SLAP documentation
    if pol
        beam_width = 18.8;
    else
        beam_width = 18.1;
    end
    footprintm = (slant*tand(beam_width/2))/cosd(elev);

    r_earth = 6378e3;
    % convert footprint radius in meters to degrees longitude
    footprint_major = footprintm / (pi*r_earth) * 180;
    % the 6-deg averaged footprint, where each individual -3 dB footprint
    % is an ellipse, is approximated as a circle
    footprint_minor = footprint_major;
    %% generate KML file to plot on Google Earth with lat&lon values and Tb/SM data
    
    % create KML file with the name specified
    fida=fopen([timestr,'_m', power, name, poltype, '.kml'],'w');
    
    line1=['<?xml version="1.0" encoding="UTF-8"?>'];
    line2=['<kml xmlns="http://earth.google.com/kml/2.0">'];
    line3=['<Document>'];
    line4=[strcat('<name>',strcat(timestr,'_m', power, name, poltype, '.kml'),'</name>')];
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
    
    % if plotting soil moisture, flip the color scale to get blues as
    % smaller values and reds as higher values
    if ~soil
        cls =     flipud(cls);
    end
    
    % to have more partitions of smaller ranges of values plotted on Google Earth
    %   cls =  ['50F00014'; '50F01414'; '50F02814'; '50F03C14'; '50F05014'; ...
    %      '50F06414'; '50F07814'; '50F08C14'; '50F0A014'; '50F0B414'; '50F0C814'; ...
    %      '50F0DC14'; '50F0F014'; '50F0FF14'; '50B4FF14'; '5078FF14'; '5014F000'; ...
    %      '5014F014'; '5014F028'; '5014F03C'; '5014F050'; '5014F064'; '5014F078'; ...
    %      '5014F08C'; '5014F0A0'; '5014F0B4'; '5014F0C8'; '5014F0DC'; '5014F0F0'; ...
    %      '5014F0FF'; '5014B4FF'; '501478FF'; '50143CFF'; '501400FF'; '501400FF'; ...
    %      '501400DC'; '501400C8'; '501400B4'; '501400A0'; '5014008C'; '50140078';'50140064'];
    
    fprintf(fida,'%s\r%s\r%s\r%s\r%s\r%s\r%s\r%s\r', ...
        line1,line2,line3,line4,line6,line7,line8,line9);
    
    %-----------------------------------
    % bin into value ranges (20 colors and 2 outside of range)
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
    % to have more partitions of smaller ranges of values plotted on Google Earth
    % f24=find(Tb>=minBin+deltaB*23 & Tb<=minBin+deltaB*24);
    % f25=find(Tb>=minBin+deltaB*24 & Tb<=minBin+deltaB*25);
    % f26=find(Tb>=minBin+deltaB*25 & Tb<=minBin+deltaB*26);
    % f27=find(Tb>=minBin+deltaB*26 & Tb<=minBin+deltaB*27);
    % f28=find(Tb>=minBin+deltaB*27 & Tb<=minBin+deltaB*28);
    % f29=find(Tb>=minBin+deltaB*28 & Tb<=minBin+deltaB*29);
    % f30=find(Tb>=minBin+deltaB*29 & Tb<=minBin+deltaB*30);
    % f31=find(Tb>=minBin+deltaB*30 & Tb<=minBin+deltaB*31);
    % f32=find(Tb>=minBin+deltaB*31 & Tb<=minBin+deltaB*32);
    % f33=find(Tb>=minBin+deltaB*32 & Tb<=minBin+deltaB*33);
    % f34=find(Tb>=minBin+deltaB*33 & Tb<=minBin+deltaB*34);
    % f35=find(Tb>=minBin+deltaB*34 & Tb<=minBin+deltaB*35);
    % f36=find(Tb>=minBin+deltaB*35 & Tb<=minBin+deltaB*36);
    % f37=find(Tb>=minBin+deltaB*36 & Tb<=minBin+deltaB*37);
    % f38=find(Tb>=minBin+deltaB*37 & Tb<=minBin+deltaB*38);
    % f39=find(Tb>=minBin+deltaB*38 & Tb<=minBin+deltaB*39);
    % f40=find(Tb>=minBin+deltaB*39 & Tb<=minBin+deltaB*40);
    % f41=find(Tb>=minBin+deltaB*40 & Tb<=minBin+deltaB*41);
    % f42=find(Tb>=minBin+deltaB*41 & Tb<=minBin+deltaB*42);
    % f43=find(Tb>maxBin);
    %
    
    
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