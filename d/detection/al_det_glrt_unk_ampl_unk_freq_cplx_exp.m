function [Tx,thresh] = al_det_glrt_unk_ampl_unk_freq_cplx_exp(x,var_wgn,Pfa)
% [Tx,thresh] = al_det_glrt_unk_ampl_unk_freq_cplx_exp(x,var_wgn,Pfa)
% GLRT Detector. It detects a complex exponent signal 
% with unknown complex amplitude and frequency.
%
% Input paramaters:
% x - input signal x,
% var_wgn - Noise Variance,
% Pfa - required Probability of the False Alarm.
%
% Output paramaters:
% Tx - Test Statistic,
% thresh - Threshold gamma for the given Pfa.
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development,
% Algorithm 12.21 – Unknown complex amplitude and frequency (also Algorithm 10.8)

% % Detector Statistic from periodogram made by FFT
N = length(x);
X = fft(x);
Xa = abs(X);
Xa2 = Xa.^2*(1/N);
Tx = max(Xa2);
% Threshold gamma
thresh = (var_wgn)*log((N-1)/Pfa);

return
