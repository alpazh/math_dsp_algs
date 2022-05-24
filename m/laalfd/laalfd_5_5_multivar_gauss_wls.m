close all
clear
clc
% fix the random seed to reproduce results
rng(12345)

%% Generate Gaussian random data for given covariance

% define number of dimensions M of vector to generate
M = 10
N = 5
% Variant 1: define unity covariance matrix 
% (in this case LS estimate is equal to WLS estimate)
% (BLUE theorem of Gauss)
% CM = eye(M);

% Variant 2: define non-diagonal covariance matrix
CM = randn(M,M)*diag(rand(1,M));
CM = CM'*CM;

% define number of data points

% generate N samples of M-dimensional random vectors
% x = randn(M,N)*1.0;
% make Cholesky factorization of the covariance matrix
F = chol(CM)

% compute correlated random vectors
eu = randn(M,1);
F = F.'
e = F*eu;

A = randn(M,N)
x = round(rand(N,1)*10)/10
b = A*x
bn = b + e

%% Compare LS to Weighted LS
x_hat_ls = inv(A'*A)*A'*bn;
x_hat_wls = inv(A'*inv(CM)*A)*A'*inv(CM)*bn;
[x x_hat_ls x_hat_wls]
e_ls  = x_hat_ls -x 
e_wls = x_hat_wls-x
norm_e_ls = norm(e_ls) 
norm_e_wls = norm(e_wls)

% covariance W of the estimate x_hat_wls
W = inv(A'*inv(CM)*A)

return
