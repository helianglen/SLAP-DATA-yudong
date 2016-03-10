%% check for 1.09 Hz noise in FFT of V-pol Data
% V2_I = squeeze(FullMomGroup.m2_ant(:,3,:));
% V2_I = reshape(V2_I', numel(V2_I),1);

v2 = v2ant;
Fs = 1/5e-4;
v2fft = fft(v2, 2^nextpow2(length(v2))) /numel(v2);
f = Fs/2*linspace(0,1,numel(v2fft)/2+1)';
figure; plot(f, 2*abs(v2fft(1:numel(v2fft)/2+1)))
ylabel('|Y(f)|')
xlabel('Frequency (Hz)')
title('Single-Sided Amplitude Spectrum of y(t)')
ylim([0 10])
xlim([0 100])