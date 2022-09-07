function [Tx,thresh] = al_det_energy_detector_cplx(x,var_wgn,Pfa)
% [Tx,thresh] = al_det_energy_detector_cplx(x,Cs,var_wgn)
% Complex Energy Detector for completely unknown signals
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
% Algorithm 12.18 – Complex Energy Detector (also Algorithm 10.5)

N = length(x);

% Detector Statistic
Tx = real(x'*x);
% Threshold gamma
x_max = N + 10*sqrt(2*N);
thresh = al_compute_energy_detector_thresh(2*N,Pfa,var_wgn/2,x_max);

return
