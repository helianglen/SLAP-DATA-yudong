
% Input: 
%	x: time series
%       deltaT: sampling interval in seconds 
% Output: 
%	f: frequency vector, in Hz 
%  	P1: Single-Sided Amplitude Spectrum 

function [f, P1] = fft_spectra (x, deltaT) 

Fs = 1 / deltaT;    %sampling freq
vlen = floor(length(x)/2)*2   % make length even
vfft = fft(x(1:vlen));
P2 = abs(vfft/vlen);
P1 = P2(1:vlen/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(vlen/2))/vlen;

end 


