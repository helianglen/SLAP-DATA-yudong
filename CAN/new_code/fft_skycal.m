%% Visualize, histogram, and fft sky calbration day data 
%  Use predefined ranges to identify foam or sky targets 
%  These ranges are a priori expert knowledge and needs to be adjusted case by case
%  We use a two-round method to get rid of outliers: 
%   The first round computes mean and std from the data
%   The 2nd round gets rid of outliers (out of 3-sigma range) and recompute mean and std 

clear

% sampling interval: 0.5 ms
 deltaT = 0.5e-3; 

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
    %load(filesRad(i).name(1:end-11), 'time')     % RADTELEM_20151005T145000_FB.mat
    load(filesRad(i).name)                       % RADTELEM_20151005T145000_FB_m2data.mat
 
    timeAll = [timeAll; time_file]; 
    h2antAll = [h2antAll; h2ant]; 
    v2antAll = [v2antAll; v2ant]; 

end


if 0 
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

end 

% there is a gap near the end, so cut it off
h2sky_ind = find ( h2antAll > skyLow & h2antAll < skyHigh ); 
h2sky=h2antAll(h2sky_ind);  
h2sky_time=timeAll(h2sky_ind); 

v2sky_ind= find ( v2antAll > skyLow & v2antAll < skyHigh ); 
v2sky=v2antAll( v2sky_ind ); 
v2sky_time=timeAll(v2sky_ind) ; 

h2foam_ind = find(h2antAll > foamLow & h2antAll < foamHigh); 
h2foam=h2antAll(h2foam_ind); 
h2foam_time=timeAll(h2foam_ind); 

v2foam_ind = find(v2antAll > foamLow & v2antAll < foamHigh); 
v2foam=v2antAll(v2foam_ind); 
v2foam_time=timeAll(v2foam_ind); 

% -------- sky data ----------------------------- 

figure
subplot(3, 2, 1)
%plot(v2sky_time, v2sky, 'b')
plot(v2sky, 'b')
title('Time series of V-pol Sky Cal, Oct. 5, 2015')
axis([0 inf 2.2e6  2.5e6]); 
%datetick('x', 15)

subplot(3, 2, 2)
%plot(h2sky_time, h2sky, 'g')
plot(h2sky, 'g')
%plot(h2sky_time, h2sky, 'g')
axis([0 inf 2.4e6  2.7e6]); 
title('Time series of H-pol Sky Cal, Oct. 5, 2015')
%datetick('x', 15)

subplot(3, 2, 3)
h1=histogram(v2sky, 2000);
axis([2.2e6  2.5e6  0 inf]); 
h1.FaceColor='b';
h1.EdgeColor='b';
title('Histograms of V-pol Sky Cal, Oct. 5, 2015')

subplot(3, 2, 4)
h2=histogram(h2sky, 2000);
axis([2.4e6  2.7e6  0 inf]); 
h2.FaceColor='g';
h2.EdgeColor='g';
title('Histograms of H-pol Sky Cal, Oct. 5, 2015');

subplot(3, 2, 5)
[vf, vP1] = fft_spectra(v2sky, deltaT); 
semilogx(vf(2:end),vP1(2:end))
%axis([0 200 0 0.5e3]); 
axis([0.1 1000 0 0.05]); 
title('Noise Power Spectrum, V-pol')
set(gca, 'XTick', [0.1 1 10 100 1000])
xlabel('f (Hz)')
ylabel('Power')

subplot(3, 2, 6)
[hf, hP1] = fft_spectra(h2sky, deltaT);
semilogx(hf(2:end),hP1(2:end))
%axis([0 200 0 0.5e3]); 
axis([0.1 1000 0 0.05]); 
title('Noise Power Spectrum, H-pol')
set(gca, 'XTick', [0.1 1 10 100 1000])
xlabel('f (Hz)')
ylabel('Power')

print('sky_ts_hist_fft.png', '-dpng');

% -------- foambox data ----------------------------- 

figure
subplot(3, 2, 1)
%plot(v2foam_time, v2foam, 'b')
plot(v2foam, 'b')
title('Time series of V-pol Foam Cal, Oct. 5, 2015')
axis([0 inf 3.2e6  3.6e6]); 
%datetick('x', 15)

subplot(3, 2, 2)
%plot(h2foam_time, h2foam, 'g')
plot(h2foam, 'g')
%plot(h2foam_time, h2foam, 'g')
axis([0 inf 3.4e6  3.8e6]); 
title('Time series of H-pol Foam Cal, Oct. 5, 2015')
%datetick('x', 15)

subplot(3, 2, 3)
h1=histogram(v2foam, 2000);
axis([3.2e6  3.6e6  0 inf]); 
h1.FaceColor='b';
h1.EdgeColor='b';
title('Histograms of V-pol Foam Cal, Oct. 5, 2015')

subplot(3, 2, 4)
h2=histogram(h2foam, 2000);
axis([3.4e6  3.8e6  0 inf]); 
h2.FaceColor='g';
h2.EdgeColor='g';
title('Histograms of H-pol Foam Cal, Oct. 5, 2015');

subplot(3, 2, 5)
[vf, vP1] = fft_spectra(v2foam, deltaT); 
semilogx(vf(2:end),vP1(2:end))
axis([0.1 1000 0 0.05]); 
title('Noise Power Spectrum, V-pol')
set(gca, 'XTick', [0.1 1 10 100 1000])
xlabel('f (Hz)')
ylabel('Power')

subplot(3, 2, 6)
[hf, hP1] = fft_spectra(h2foam, deltaT);
semilogx(hf(2:end),hP1(2:end))
axis([0.1 1000 0 0.05]); 
%axis([0 200 0 0.5e3]); 
title('Noise Power Spectrum, H-pol')
set(gca, 'XTick', [0.1 1 10 100 1000])
xlabel('f (Hz)')
ylabel('Power')

print('foam_ts_hist_fft.png', '-dpng');


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


