script
close all
clear global
clc
% Kay, Estimation Theory, Chapter 11, p.352, Example 11.3:
% Estimate DC level in AWGN assuming uniform prior PDF.

% rng(123)
N = 1e3;
A = 1;
A0 = A;
sigma = 1;
mu_A = 0;
sigma_A = A/3;
w = randn(N,1)*sigma;
x = ones(N,1)*A + w;

A_hat = mean(x)
A_hat_mmse = al_est_bayes_mmse_dc_level_in_wgn(x,sigma,sigma_A,mu_A)
A_hat_map = al_est_bayes_map_dc_level_in_wgn_uni_prior(x,A0)

err = A_hat - A
err_mmse = A_hat_mmse - A
err_map = A_hat_map - A

return
