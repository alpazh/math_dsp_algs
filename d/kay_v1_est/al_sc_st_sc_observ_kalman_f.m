function [s_hat_n_n_v,M_n_n_v,M_n_n_1_v,MMSE_v] = ...
    al_sc_st_sc_observ_kalman_f(x,sigma_n,sigma_u,a,M_n_1_n_1,s_hat_n_1_n_1)
% [s_hat_n_n_v,M_n_n_v,M_n_n_1_v,MMSE_v] = al_sc_st_sc_observ_kalman_f(x,sigma_n,sigma_u,a,M_n_1_n_1,s_hat_n_1_n_1)  
% Scalar state scalar observation Kalman Filter
% Input arguments:
% x - observed signal
% sigma_n - variance of observation noise (time-varying)
% sigma_u - - variance of excitation noise
% a - Gauss Markov model parameter
% M_n_1_n_1 - initial MSE
% s_hat_n_1_n_1 - initial estimate
% Output arguments:
% s_hat_n_n_v - Estimated signal (Kalman Filter Output)
% M_n_n_v - Kalman Filter MSE
% M_n_n_1_v - Prediction MSE
% MMSE_v - vector of interleaved Kalman Filter MSE and Prediction MSE
% Reference: Kay, Estimation Theory, Chapter 12.6 Sequential LMMSE
% Estimation.

%   Author: Alexey Zherebtsov

N = length(x);
s_hat_n_n_v = zeros(N,1);
M_n_n_v = zeros(N,1);
M_n_n_1_v = zeros(N,1);
MMSE_v = zeros(2*N,1);
MMSE_v(2*1) = M_n_1_n_1;
M_n_n_v(1) = M_n_1_n_1;
for k = 2:N
    % Prediction
    s_hat_n_n_1 = a*s_hat_n_1_n_1;
    % Min Prediction MSE
    M_n_n_1 = a^2*M_n_1_n_1 + sigma_u^2;
    M_n_n_1_v(k) = M_n_n_1;
    MMSE_v(2*k-1) = M_n_n_1;
    % Kalman Gain
    K_n = M_n_n_1/(sigma_n(k)^2 + M_n_n_1);
    % Correction
    s_hat_n_n = s_hat_n_n_1 + K_n * (x(k) - s_hat_n_n_1);
    s_hat_n_1_n_1 = s_hat_n_n;
    % Min MSE
    M_n_n = (1 - K_n)*M_n_n_1;
    M_n_1_n_1 = M_n_n;
    s_hat_n_n_v(k) = s_hat_n_n;
    M_n_n_v(k) = M_n_n;
    MMSE_v(2*k) = M_n_n;
end

end
