%%
% uses sum and mean insteald of nansum and nanmean to avoid license restrictions
clear

%% Initialize key options
flagHalf = 0;
% flag to process data for fore or aft half-scan
fore = 1;

% flag to plot h or v pol data on Google Earth
GoogleEarth  = 1;

%% roll flag and limits
flagRoll = 1;
% Aircraft roll angle limits if removing data during a turn
rollmin = -5;
rollmax = 5;

% RDRRETURN_20151030T181002.mat
files = dir('6deg_RDRRETURN_*.mat');
files_to_process = 1:length(files);
%%
% beam width information by polarization from SLAP documentation
beam_width_h_pol = 20.5;
beam_width_v_pol = 19.5;

[time_total, az_total, lat_total, lon_total, hh_total, vv_total, hv_total, ...
 vh_total, alt_total, trk_total, roll_total, hdg_total] = deal([]);


%%
for k = files_to_process
    files(k).name
    % load time, radar_h, and radar_v data
    load(files(k).name)
    
    time_total = [time_total; timeavg];
    hh_total = [hh_total; mean_hh];
    vv_total = [vv_total; mean_vv];
    hv_total = [hv_total; mean_hv];
    vh_total = [vh_total; mean_vh];
    az_total = [az_total; azavg];
    lat_total = [lat_total; latavg];
    lon_total = [lon_total; lonavg];
    alt_total = [alt_total; altavg];
    trk_total = [trk_total; trkavg];
    roll_total = [roll_total; rollavg];
    hdg_total = [hdg_total; hdgavg];

end   % file loop 
    %-------------------------------------------------------------------
    if GoogleEarth


       % only plot on either the fore or aft half scan if flagHalf is on      
        if flagHalf
            % initialize flag_half_scan as zeros
            flag_half_scan = zeros(length(trk_total),1);
            
            for k = 1:length(trk_total)
                % get 180 degrees of scan width around track angle to get front half
                scan_lower = trk_total(k) - 90;
                scan_upper = trk_total(k) + 90;
                scan_upper(scan_upper>360) = scan_upper(scan_upper>360) - 360;
                scan_lower(scan_lower<0) = scan_lower(scan_lower<0) + 360;
                
                % get indices of azimuth angles within the 180 degrees of the front
                % half
                flagf1 = az_total(k)<scan_upper;
                flagf2 = az_total(k)>=scan_lower;
                
                if trk_total(k) < 180;
                    if scan_lower > 180 && (flagf1 || flagf2)
                        flag_half_scan(k) = k;
                    elseif scan_lower <= 180 && (flagf1 && flagf2)
                        flag_half_scan(k) = k;
                    end
                else
                    if scan_upper > 180 && (flagf1 && flagf2)
                        flag_half_scan(k) = k;
                    elseif scan_upper <= 180 && (flagf1 || flagf2)
                        flag_half_scan(k) = k;
                    end
                end
            end
            
            % turn flag_half_scan into a vector of logicals where ones are indices of
            % fore half scan data points
            flag_half_scan = flag_half_scan>0;
            
            if ~fore
                % get indices of aft half scan as the opposite of the indices of the fore
                % half scan if the flag for the fore half scan is turned off
                flag_half_scan = ~flag_half_scan;
            end
            
        else
            % if flagHalf is off, to plot fore and aft half scans together, get all logical ones
            flag_half_scan = find(lat_total>0);
        end
        
        % get geolocation data for fore or aft half scans using flag_half_scan
        trk_half_scan = trk_total(flag_half_scan);
        lon_half_scan = lon_total(flag_half_scan);
        lat_half_scan = lat_total(flag_half_scan);
        roll_half_scan = roll_total(flag_half_scan);
        vv_half_scan = vv_total(flag_half_scan);
        hh_half_scan = hh_total(flag_half_scan);
        hv_half_scan = hv_total(flag_half_scan);
        vh_half_scan = vh_total(flag_half_scan);
        alt_half_scan = alt_total(flag_half_scan);
        time_half_scan = time_total(flag_half_scan);
        az_half_scan = az_total(flag_half_scan);
        hdg_half_scan = hdg_total(flag_half_scan);
        
        Tb_hh = hh_half_scan;
        Tb_hv = hv_half_scan;
        Tb_vv = vv_half_scan;
        Tb_vh = vh_half_scan;
        
        Lat = lat_half_scan;
        Long = lon_half_scan;
        Alt = alt_half_scan;
        Roll = roll_half_scan; 
        

        % if flagRoll is set on, perform the roll filter
        if flagRoll
            indRoll = find(roll_half_scan < rollmax & roll_half_scan > rollmin);
            Tb_hv = Tb_hv(indRoll);
            Tb_vh = Tb_vh(indRoll);
            Tb_vv = Tb_vv(indRoll);
            Tb_hh = Tb_hh(indRoll);
            Lat = lat_half_scan(indRoll);
            Long = lon_half_scan(indRoll);
            Alt = alt_half_scan(indRoll);
            Roll = roll_half_scan(indRoll); 
        end

        minBin= - 15;
        maxBin=  15;
        
        'about to plot ...'
        %keyboard
        %save_kml('radar_Fore_HV.kml', Tb_hv, Lat, Long, Alt, Roll, beam_width_h_pol, minBin, maxBin) 
        %save_kml('radar_Fore_VH.kml', Tb_vh, Lat, Long, Alt, Roll, beam_width_v_pol, minBin, maxBin) 
        save_kml('radar_Fore_HH.kml', Tb_hh, Lat, Long, Alt, Roll, beam_width_h_pol, minBin, maxBin) 
        save_kml('radar_Fore_VV.kml', Tb_vv, Lat, Long, Alt, Roll, beam_width_v_pol, minBin, maxBin) 

        
    end  % if Google Earth 
    %-------------------------------------------------------------------
quit

