function [Tx,s_hat] = al_det_estimator_correlator_cplx(x,Cs,var_wgn)
% [Tx,s_hat] = al_det_estimator_correlator_cplx(x,Cs,var_wgn)
% Complex Estimator Correlator Detector for deterministic signals
% Threshold is hard to find analytically and is omitted.
%
% Input paramaters:
% x - input signal x,
% Cs - Random Signal Covariance Matrix,
% var_wgn - Noise Variance.
%
% Output paramaters:
% Tx - Test Statistic,
% s_hat - estimated random signal s.
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development,
% Algorithm 12.17 – Complex Estimator-Correlator Detector (also Algorithm 10.4)

N = size(Cs,1);
% Estimated Signal
s_hat = Cs*inv(Cs+var_wgn*eye(N))*x;
% Detector Statistic
Tx = real(x'*s_hat);

return
