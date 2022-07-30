script
close all
clear global
clc
% rng(123)
N = 1e3;
A = 1;
sigma = 1;
w = randn(N,1)*sigma;
x = ones(N,1)*A + w;

A_hat = mean(x);
err = A_hat - A;

return
