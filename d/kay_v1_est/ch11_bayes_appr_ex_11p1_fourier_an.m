script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 11, p.347, Example 11.1:
% Estimate amplitude of sum sin+cos signals of the same frequency in the
% AWGN.

% rng(123)
N = 1e3;
sigma = 9;
a = 1;
b = 3;
Fs = 1
f0 = 0.05*Fs
t = (0:1/Fs:(N-1)/Fs)';
w = randn(N,1)*sigma;
x = a*cos(2*pi*f0*t) + b*sin(2*pi*f0*t) + w;

figure
plot(x,'b- .'),grid on,hold on

sigma = 1;
sigma_theta = 1;
theta_hat = al_est_bayes_mmse_sin_cos_magn_est(x,sigma,sigma_theta,f0,Fs)

if N < 1e4
    mu_theta = [0 0]';
    C_theta = eye(2)*sigma_theta^2;
    C_w = eye(length(w))*(sigma^2);
    H = [cos(2*pi*f0*t) sin(2*pi*f0*t)];
    size(C_w)
    theta_hat_gen_eq = al_est_bayes_mmse_linear_model(x,H,mu_theta,C_theta,C_w)
end
theta_hat
return
