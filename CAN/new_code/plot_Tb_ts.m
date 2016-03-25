%% Visualize calbrated Tb data 
% 
%   The 2nd round gets rid of outliers (out of 3-sigma range) and recompute mean and std 

clear

filesRad = dir('RAD*FB_m2data.mat');
files_to_process = 1:length(filesRad);

timeAll = [];
h2antAll = [];
v2antAll = [];

for i = files_to_process

    % time data from FB.mat and h2ant and v2ant from FB_m2data.mat 
    load(filesRad(i).name)                       % RADTELEM_20151005T145000_FB_m2data.mat
 
    timeAll = [timeAll; time_file]; 
    h2antAll = [h2antAll; h2ant]; 
    v2antAll = [v2antAll; v2ant]; 

end

       % Two-point calibration fit form Sky Cal day data
        % it seems not working: for 20151030T171645, got Tbh > Tbv
        % countH = 2507448.51 + TbH * 4152.4488
        % countV = 2306412.72 + TbV * 4376.8796
        tbcalH = (h2antAll - 2507448.51) / 4152.4488;
        tbcalV = (v2antAll - 2306412.72) / 4376.8796;

figure
subplot(2, 1, 1) 
plot(timeAll, tbcalH, 'g')
axis([min(timeAll) max(timeAll) 50 450]) 
title('H (green) and V (blue) Tb (K)')
datetick('x', 15)
hold on
s=plot(timeAll, tbcalV, 'b')

subplot(2, 1, 2) 
plot(timeAll, tbcalV-tbcalH, 'Color', [0.5 0.5 0.5])
axis([min(timeAll) max(timeAll) -200 200]) 
title('TbV - TbH (K)')
datetick('x', 15)

print('calTb-ts.png', '-dpng');  

