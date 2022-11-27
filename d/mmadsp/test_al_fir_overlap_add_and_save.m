%script
close all
clear
clc
% Test for the fast FFT-based filters
% author: Aleksei Zherebtsov

if(exist('OCTAVE_VERSION', 'builtin'))
    pkg load communications
    pkg load signal
end


% sampling frequency
Fs=8e3; 
% signal length
Np=3e3; 
% first sine signal frequency
F1=233; 
% second sine signal frequency
F2=415; 
% third sine signal frequency
F3=1877; 
% sampling period
Ts = 1/Fs; 
% time vector
t = (Ts:Ts:Np*Ts).'; 
% generate input signal
s = sin(2*pi*F1*t) + sin(2*pi*F2*t) + sin(2*pi*F3*t);
% add noise
s = s + wgn(Np,1,0.01);
% filter order
FiltOrd = 384;
% FFT size
Nfft = 1024;

% design filter and save filter coefficients
h = fir1(FiltOrd,1/4).';

% filter signal in time-domain to check the correctness of the overlap-add algorithm
s_td = filter(h,1,s);
% filter signal in frequency-domain by the overlap-add algorithm
s_oa = al_fir_overlap_add(h,s,Nfft);
% filter signal in frequency-domain by the overlap-save algorithm
s_os = al_fir_overlap_save(h,s,Nfft);

% check that all three algorithms provide the same results
figure
plot(s_td, 'b- s'),grid on,hold on
plot(s_oa, 'r- o'),grid on,hold on
plot(s_os, 'g- x'),grid on,hold off
legend({'Time-Domain','Overlap-Add','Overlap-Save'})

% check that the filtered signal spectrum does not contain high frequency components
figure
[PSDPSD,f]=pwelch(s_oa,[],[],[],Fs,'twosided');
plot(f-Fs/2,fftshift(10*log10(PSDPSD/max(PSDPSD))),'b- '),grid on,title('Spectrum of the signal filtered by the Overlap-Add algorithm');

return
