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
    load(filesRad(i).name(1:end-11), 'time')     % RADTELEM_20151005T145000_FB.mat
    load(filesRad(i).name)                       % RADTELEM_20151005T145000_FB_m2data.mat
 
    timeAll = [timeAll; time]; 
    h2antAll = [h2antAll; h2ant]; 
    v2antAll = [v2antAll; v2ant]; 

end

figure
subplot(2, 2, 1) 
plot(v2antAll, 'b')
title('Time series of V-pol counts, Oct. 5, 2015')

subplot(2, 2, 2) 
plot(h2antAll, 'g')
title('Time series of H-pol counts, Oct. 5, 2015')

subplot(2, 2, 3) 
h1=histogram(v2antAll, 2000);
h1.FaceColor='b';
h1.EdgeColor='b';
title('Histograms of V-pol counts, Oct. 5, 2015')


subplot(2, 2, 4) 
h2=histogram(h2antAll, 2000); 
h2.FaceColor='g';
h2.EdgeColor='g';
title('Histograms of H-pol counts, Oct. 5, 2015'); 

print('ts_hist.png', '-dpng');  

h2sky=h2antAll(h2antAll > skyLow & h2antAll < skyHigh); 
v2sky=v2antAll(v2antAll > skyLow & v2antAll < skyHigh); 

h2foam=h2antAll(h2antAll > foamLow & h2antAll < foamHigh); 
v2foam=v2antAll(v2antAll > foamLow & v2antAll < foamHigh); 

disp('mean(v2foam) mean(h2foam) mean(v2sky) mean(h2sky) std(v2foam) std(h2foam) std(v2sky) std(h2sky)') 
disp('----------------- before outlier removal ---------------------------------------------------------------')
fprintf('%10.0f, ', [mean(v2foam), mean(h2foam), mean(v2sky), mean(h2sky), ...
       std(v2foam),  std(h2foam),  std(v2sky),  std(h2sky)]); 
fprintf('\n');  

%new ranges 
new_h2skyLow = mean(h2sky) - 3.0 * std(h2sky);  
new_h2skyHigh = mean(h2sky) + 3.0 * std(h2sky);  

new_v2skyLow = mean(v2sky) - 3.0 * std(v2sky); 
new_v2skyHigh = mean(v2sky) + 3.0 * std(v2sky); 

new_h2foamLow = mean(h2foam) - 3.0 * std(h2foam);
new_h2foamHigh = mean(h2foam) + 3.0 * std(h2foam);  

new_v2foamLow = mean(v2foam) - 3.0 * std(v2foam);
new_v2foamHigh = mean(v2foam) + 3.0 * std(v2foam);

new_h2sky=h2antAll(h2antAll > new_h2skyLow & h2antAll < new_h2skyHigh); 
new_v2sky=v2antAll(v2antAll > new_v2skyLow & v2antAll < new_v2skyHigh); 

new_h2foam=h2antAll(h2antAll > new_h2foamLow & h2antAll < new_h2foamHigh); 
new_v2foam=v2antAll(v2antAll > new_v2foamLow & v2antAll < new_v2foamHigh); 

disp('----------------- after outlier removal ---------------------------------------------------------------')
fprintf('%10.0f, ', [mean(new_v2foam), mean(new_h2foam), mean(new_v2sky), mean(new_h2sky), ...
       std(new_v2foam),  std(new_h2foam),  std(new_v2sky),  std(new_h2sky)]); 
fprintf('\n');  


