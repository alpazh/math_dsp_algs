function [Tx,thresh] = al_det_replica_correlator_cplx(s,x,var_wgn,Pfa)
% [Tx,thresh] = al_det_replica_correlator_cplx(s,x,var_wgn,Pfa)
% Complex Replica-Correlator Detector for deterministic signals
%
% Input paramaters:
% s - known deterministic signal to be detected,
% x - input signal x,
% var_wgn - known WGN variance
% Pfa - required Probability of the False Alarm
%
% Output paramaters:
% Tx - Test Statistic,
% thresh - THreshold gamma for the given Pfa.
% Test example of usage is in d\detection\test_det_ml_pam_odd.m
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development,
% Algorithm 12.14 – Replica-correlator (matched filter) (also Algorithm 10.1)


% Detector Statistic
Tx = real(x'*s);
% gamma threshold
thresh = sqrt((var_wgn/2)*(s'*s))*al_q_inv_func(Pfa);

return
