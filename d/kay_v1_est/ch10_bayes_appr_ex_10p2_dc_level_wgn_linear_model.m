script
close all
clear global
clc
% rng(123)
N = 3e0;
A = 1.234;
sigma = 1;
mu_A = 0;
sigma_A = A/3;
w = randn(N,1)*sigma;
x = ones(N,1)*A + w;

A_hat = mean(x)
A_hat_mmse = al_est_bayes_mmse_dc_level_in_wgn(x,sigma,sigma_A,mu_A)
H = ones(size(x));
mu_theta = 0;
C_theta = sigma_A^2;
C_w = eye(length(w))*(sigma^2);
size(C_w)
A_hat_mmse_lm = al_est_bayes_mmse_linear_model(x,H,mu_theta,C_theta,C_w)
A_hat_mmse

err = A_hat - A
err_mmse = A_hat_mmse - A

return
