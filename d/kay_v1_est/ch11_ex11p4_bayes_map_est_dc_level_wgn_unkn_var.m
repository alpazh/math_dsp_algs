script
close all
clear global
clc
% Kay, Estimation Theory, Chapter 11, p.355, Example 11.4:
% Estimate DC level in AWGN with unknown variance assuming gamma prior PDF.

% rng(123)
N = 1e3;
A = 3;
sigma = 2;
mu_A = 0;
gamma = 1;
alpha = 1;
w = randn(N,1)*sigma;
x = ones(N,1)*A + w;

A_hat = mean(x)
[A_hat_map,sigma_2_hat] = al_est_bayes_map_dc_level_in_wgn_uni_prior(x,gamma,alpha,mu_A)

err_A = A_hat - A
err_A_map = A_hat_map - A

err_sigma_2 = sigma_2_hat - sigma^2





return
