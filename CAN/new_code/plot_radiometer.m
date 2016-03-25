clear
% flag to choose whether to perform calculations to put data in tabular form
data_accum = 0;
% flag to plot h or v pol data on Google Earth
GoogleEarth  = 1;
pol = 1; % 1 == hpol and 0 == vpol
% flag to process data for fore or aft half-scan
fore = 1;
% flag to decide whether to plot all of the scan or just half, specifically
% the fore or aft halfs
flagHalf = 0;
% load RAD files 
filesRad2 = dir('RAD*FB_v4.mat');
%YDT filesRad2 = dir('RADTELEM_20151108T203000_FB_v4.mat');


% choose whether to plot only NE or SW facing part of flight lines for May
% 2014 iPHEX Campaign AM flights over Western NC. If the NEorSW variable is enabled,
% the plotting algorithm chooses one side or another to plot based on flag (NE = 1, SW = 0).
% NEorSW = 1;

% flag to remove streaks in h-pol Tb data
removestreaks = 0;

% flag to calculate or not calculate SM
soil = 0;

%YDT, set a random date for now 
flightdatestr = pwd;
%flightdate = str2double((flightdatestr(end-8:end-7)));

flightdate = 2; 


% which type of data used in processing, currently M2
power = '2';

% Aircraft roll angle limits if removing data during a turn
flagRoll = 1;
rollmin = -4;
rollmax = 4;

% initiate variables
[SM_total, az_total, emiss_surf_total, clay_total, emissinit_total, temphswitch, tempvswitch, gamma_total, alt_total roll_total trk_total hdg_total tau_total emiss_soil_total emiss_h_total  Lat_total, ndvi_total, tbh_total, tbv_total, Lon_total, ndvi_total, lcclass_total, lst_total, time_total] = deal([]);
% save SMvals SM_total emissinit_total Lat_total gamma_total tau_total emiss_surf_total emiss_soil_total emiss_h_total Lon_total lcclass_total ndvi_total  lst_total  time_total

load daq_all

%New DAQ columns per Albert
%TempLNAH is the IMA H Amplifier, TempLNAV is the IMA V amplifier, 
%TempAntennaH is the microstrip H 9-way and TempAntennaV is the microstrip V 9-way

[~, daq_sort_i] = sort(timeDaq); 
daq_sort_i = daq_sort_i([diff(timeDaq(daq_sort_i)) ~= 0; true]);

new_timeDaq = timeDaq(daq_sort_i) 
new_IMA_H_SWITCH = IMA_H_SWITCH(daq_sort_i);
new_IMA_V_SWITCH = IMA_V_SWITCH(daq_sort_i);
new_TempDiplexerH = TempDiplexerH(daq_sort_i);
new_TempDiplexerV = TempDiplexerV(daq_sort_i);
new_TempLNAH = TempLNAH(daq_sort_i);
new_TempLNAV = TempLNAV(daq_sort_i);
new_TempAntennaH = TempAntennaH(daq_sort_i);
new_TempAntennaV = TempAntennaV(daq_sort_i); 


% =======================================================
for i = 1:length(filesRad2)
    i
    clear noscan
        load(filesRad2(i).name)
%         azavg = azavg + hdgavg;
        filesRad2(i).name
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
        'inside 0 ....'
        for k = 1:length(trkavg)
            'inside 1 ....'
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

        
        % Two-point calibration fit form Sky Cal day data 
        % it seems not working: for 20151030T171645, got Tbh > Tbv 
        % countH = 2507448.51 + TbH * 4152.4488
        % countV = 2306412.72 + TbV * 4376.8796
        tbcalh = (h_half_scan - 2507448.51) / 4152.4488; 
        tbcalv = (v_half_scan - 2306412.72) / 4376.8796; 

        % old fashion
        % skycalh = 2548245
        % skycalv = 2341245
        % tb_skycal = 10  %K

        % counts_foambox_h = 3671527
        % counts_foambox_v = 3506909
        % tb_foambox = 280.8

        % tbcalh = interp1([skycalh counts_foambox_h], [tb_skycal,  tb_foambox], h_half_scan, 'linear', 'extrap');
        % tbcalv = interp1([skycalv counts_foambox_v], [tb_skycal,  tb_foambox], v_half_scan, 'linear', 'extrap');

        'check calibration '
        %keyboard
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
        
        'here 1' 
        %keyboard
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

        'here 1.5' 
        %keyboard
        
        temp_dip_interp    = interp1(new_timeDaq,new_TempDiplexerH, time_half_scan);
        temp_lna_interp    = interp1(new_timeDaq,new_TempLNAH, time_half_scan);
        temp_ant_interp    = interp1(new_timeDaq,new_TempAntennaH, time_half_scan);
        temp_switch_interp = interp1(new_timeDaq,new_IMA_H_SWITCH, time_half_scan);
        
        temp_dip_interp(isnan(temp_dip_interp)) = nanmean(temp_dip_interp);
        temp_lna_interp(isnan(temp_lna_interp)) = nanmean(temp_lna_interp);
        temp_ant_interp(isnan(temp_ant_interp)) = nanmean(temp_ant_interp);
        temp_switch_interp(isnan(temp_switch_interp)) = nanmean(temp_switch_interp);
        temp_switch_interph = temp_switch_interp;
        
        Tb = Tb*(1- loss1) + (loss1)*( temp_ant_interp+273.16);
        Tb = Tb*(1- loss2) + (loss2)*( temp_dip_interp+273.16);
        Tb = Tb*(1- loss3) + (loss3)*( temp_dip_interp+273.16);
        Tbh = Tb*(1- loss4) + (loss4)*( temp_lna_interp+273.16);
        
        'check h correction' 
        %keyboard
        % do corrections for V-Pol
        Tb = tbv_half_scan;
        tb_otherpol = tbh_half_scan;
        Tb = Tb*(1 - 0.01) + tb_otherpol*0.01;
        
        temp_dip_interp    = interp1(new_timeDaq,new_TempDiplexerV, time_half_scan);
        temp_lna_interp    = interp1(new_timeDaq,new_TempLNAV, time_half_scan);
        temp_ant_interp    = interp1(new_timeDaq,new_TempAntennaV, time_half_scan);
        temp_switch_interp = interp1(new_timeDaq,new_IMA_V_SWITCH, time_half_scan);
        
        'here 2.5' 
        %keyboard

        temp_dip_interp(isnan(temp_dip_interp)) = nanmean(temp_dip_interp);
        temp_lna_interp(isnan(temp_lna_interp)) = nanmean(temp_lna_interp);
        temp_ant_interp(isnan(temp_ant_interp)) = nanmean(temp_ant_interp);
        temp_switch_interp(isnan(temp_switch_interp)) = nanmean(temp_switch_interp);
        temp_switch_interpv = temp_switch_interp;
        
        Tb = Tb*(1- loss1) + (loss1)*( temp_ant_interp+273.16);
        Tb = Tb*(1- loss2) + (loss2)*( temp_dip_interp+273.16);
        Tb = Tb*(1- loss3) + (loss3)*( temp_dip_interp+273.16);
        Tbv = Tb*(1- loss4) + (loss4)*( temp_lna_interp+273.16);
        
        'check v correction ' 
        %keyboard
        
        % accumulate all 10 minute data sets into one vector for each variable
        tbh_total = [tbh_total; Tbh];
        tbv_total = [tbv_total; Tbv];
        %time_total = [time_total; datevec(time_half_scan)];
        time_total = [time_total; time_half_scan];
        Lat_total = [Lat_total; lat_half_scan];
        Lon_total = [Lon_total; lon_half_scan];
        alt_total = [alt_total; alt_half_scan];
        roll_total = [roll_total; roll_half_scan ];
        trk_total = [trk_total; trk_half_scan ];
        hdg_total = [hdg_total; hdg_half_scan ];
        az_total = [az_total; az_wo_hdg];
        
        temphswitch = [temphswitch; temp_switch_interph];
        tempvswitch = [tempvswitch; temp_switch_interpv];
        
        'here 4' 
        %keyboard
        
end   
% =======================================================

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


   % create a flag of all ones to plot all data
    flagdir = Lat_total > 0;
    
    Tbh = tbh_total(flagdir);
    Tbv = tbv_total(flagdir);
    time_filtered = time_total(flagdir); 
    lat_filtered = Lat_total(flagdir);
    lon_filtered = Lon_total(flagdir);
    alt_filtered = alt_total(flagdir);
    roll_filtered = roll_total(flagdir);
    hdg_filtered = hdg_total(flagdir);


    % specify what data set is current being worked on to append to filenames
    % of output
    if fore
       look_dir = 'Fore';
    else 
       look_dir = 'Aft';
    end

 if GoogleEarth

    save_kml([look_dir, '_h', '.kml'], Tbh, lat_filtered, lon_filtered, ...
                                   alt_filtered, roll_filtered, 'h') 

    save_kml([look_dir, '_v', '.kml'], Tbv, lat_filtered, lon_filtered, ...
                                   alt_filtered, roll_filtered, 'v') 

  end   % GoogleEarth


figure
subplot(3, 1, 1)
plot(time_filtered, Tbh, 'g')
axis([min(time_filtered) max(time_filtered) 150 320])
title('H (green) and V (blue) Tb (K)')
datetick('x', 15)
hold on
plot(time_filtered, Tbv, 'b')

subplot(3, 1, 2)
plot(time_filtered, Tbv-Tbh, 'Color', [0.0 0.0 0.0])
axis([min(time_filtered) max(time_filtered) -40 40])
title('TbV - TbH (K)')
datetick('x', 15)
%hold on 
%plot(time_filtered, zeros(length(time_filtered)), 'g') 

subplot(3, 1, 3)
plot(time_filtered, roll_filtered, 'Color', [0.0 0.0 0.0])
axis([min(time_filtered) max(time_filtered) -50 10])
title('Roll angle (deg)')
datetick('x', 15)

print('calTb-roll-ts.png', '-dpng');


quit
