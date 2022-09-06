function [Tx,thresh] = al_det_limiter_replica_correlator(A,x,var_n,Pfa)
% [Tx,thresh] = al_det_limiter_replica_correlator(A,x,var_wgn,Pfa)
% Limiter-Replica-Correlator Detector for deterministic signals in IID Laplacian Noise
% It detects a known signal in IID non Gaussian noise.
%
% Input paramaters:
% A - the magnitude of the known deterministic signal to be detected,
% x - input signal x,
% var_n - known variance of the Laplacian Noise,
% Pfa - required Probability of the False Alarm.
%
% Output paramaters:
% Tx - Test Statistic,
% thresh - Threshold gamma for the given Pfa.
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development,
% Algorithm 12.15 – Limiter-Replica-Correlator (also Algorithm 10.2)

N = length(x);
% Limiter output - function g(x)
g_x = sqrt(2/var_n)*sign(x);
% Detector Statistic
Tx = sum(g_x'*A);
% Signal Energy
E = N*A.^2;
% gamma threshold
thresh = sqrt((2/var_n)*E)*al_q_inv_func(Pfa);

return
