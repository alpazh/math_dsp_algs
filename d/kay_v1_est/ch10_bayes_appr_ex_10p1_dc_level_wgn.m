script
close all
clear global
clc
% rng(123)
N = 1e3;
A = 1;
sigma = 1;
mu_A = 0;
sigma_A = A/3;
w = randn(N,1)*sigma;
x = ones(N,1)*A + w;

A_hat = mean(x)
A_hat_mmse = al_est_bayes_mmse_dc_level_in_wgn(x,sigma,sigma_A,mu_A)

err = A_hat - A
err_mmse = A_hat_mmse - A

return
