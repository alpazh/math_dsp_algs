function [Tx,thresh] = al_det_glrt_kn_sig_unk_ampl_cplx(s,x,var_wgn,Pfa)
% [Tx,thresh] = al_det_glrt_kn_sig_unk_ampl_cplx(s,x,var_wgn,Pfa)
% GLRT Detector. It detects a known deterministic signal with unknown amplitude.
%
% Input paramaters:
% s - known signal to be detected,
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
% Algorithm 12.19 – Unknown complex amplitude (also Algorithm 10.6)

% Detector Statistic
xh_s = x'*s;
% xh_s*xh_s'
Tx = (xh_s*xh_s')/((var_wgn/2)*(s'*s));
% pause

% Threshold gamma
% thresh = al_q_func(Pfa/2).^2;%real valued signal case
thresh = 2*log(1/Pfa);

return
