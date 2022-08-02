script
close all
clear global
clc
% Kay, Estimation Theory, Chapter 11.5 MAP estimators, p.351, Example 11.2:
% Estimate exponential PDF parameter theta

% rng(123)
N = 1e4;
A = 1;
theta = 3;
mu = 1/theta
x = random('Exponential',mu,N,1);

figure
plot(x,'b- .'),grid on,hold on

figure
histogram(x)

% Method of moments estimator
theta_hat_mom = 1/mean(x)

% Bayesian MAP estimator
lambda = theta/3; % prior PDF(theta) parameter
theta_hat_map = 1/(mean(x)+lambda/N)
err_mom = theta_hat_mom - theta
err_map = theta_hat_map - theta
return
