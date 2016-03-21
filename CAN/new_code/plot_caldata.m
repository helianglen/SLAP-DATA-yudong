%% Visualize calbration data and compute mean/stdev 
%  Use predefined ranges to identify foam or sky targets 
%  These ranges are a priori expert knowledge and needs to be adjusted case by case
%  We use a two-round method to get rid of outliers: 
%   The first round computes mean and std from the data
%   The 2nd round gets rid of outliers (out of 3-sigma range) and recompute mean and std 

clear

% valid ranges 
skyHigh = 3e6; 
skyLow = 2e6; 
foamHigh = 4e6; 
foamLow = 3e6; 


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

figure
subplot(2, 1, 1) 
plot(timeAll, v2antAll, 'b')
axis([min(timeAll) max(timeAll) foamLow foamHigh]) 
title('Time series of V-pol counts')
datetick('x', 15)

subplot(2, 1, 2) 
plot(timeAll, h2antAll, 'g')
axis([min(timeAll) max(timeAll) foamLow foamHigh]) 
title('Time series of H-pol counts')
datetick('x', 15)

print('ts.png', '-dpng');  

h2foam=h2antAll(h2antAll > foamLow & h2antAll < foamHigh); 
v2foam=v2antAll(v2antAll > foamLow & v2antAll < foamHigh); 

disp('mean(v2foam) mean(h2foam) mean(v2sky) mean(h2sky) std(v2foam) std(h2foam) std(v2sky) std(h2sky)') 
disp('----------------- before outlier removal ---------------------------------------------------------------')
fprintf('%10.0f, ', [mean(v2foam), mean(h2foam), NaN, NaN, ...
       std(v2foam),  std(h2foam),  NaN,  NaN]); 
fprintf('\n');  

%new ranges 

new_h2foamLow = mean(h2foam) - 3.0 * std(h2foam);
new_h2foamHigh = mean(h2foam) + 3.0 * std(h2foam);  

new_v2foamLow = mean(v2foam) - 3.0 * std(v2foam);
new_v2foamHigh = mean(v2foam) + 3.0 * std(v2foam);

new_h2foam=h2antAll(h2antAll > new_h2foamLow & h2antAll < new_h2foamHigh); 
new_v2foam=v2antAll(v2antAll > new_v2foamLow & v2antAll < new_v2foamHigh); 

disp('----------------- after outlier removal ---------------------------------------------------------------')
fprintf('%10.0f, ', [mean(new_v2foam), mean(new_h2foam), NaN, NaN, ...
       std(new_v2foam),  std(new_h2foam),  NaN,  NaN]); 
fprintf('\n');  


