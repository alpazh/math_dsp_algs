function [Tx,thresh] = al_det_glrt_unk_ampl_cplx_exp(f0,x,var_wgn,Pfa)
% [Tx,thresh] = al_det_glrt_unk_ampl_cplx_exp(f0,x,var_wgn,Pfa)
% GLRT Detector. It detects complex exponent with known frequency and 
% unknown amplitude.
%
% Input paramaters:
% f0 - known frequency of the complex signal to be detected,
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
% Algorithm 12.20 – Unknown complex amplitude of complex exponent (also Algorithm 10.7)

% Detector Statistic
nv = (0:length(x)-1)';
s = exp(1j*2*pi*f0*nv);
xh_s = x'*s;
% xh_s*xh_s'
Tx = (xh_s*xh_s')/(var_wgn/2);
% pause

% Threshold gamma
% thresh = al_q_func(Pfa/2).^2;%real valued signal case
thresh = 2*log(1/Pfa);

return
