script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 12, p.399, Example 12.3:
% Estimate amplitude of sum sin+cos signals of the same frequency in the
% AWGN.

% rng(123)
N = 5e2;
sigma = 3;
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
theta_hat_bayes_mmse = al_est_bayes_mmse_sin_cos_magn_est(x,sigma,sigma_theta,f0,Fs)
if N < 1e4
    mu_theta = [0 0]';
    C_theta = eye(2)*sigma_theta^2;
    C_w = eye(length(w))*(sigma^2);
    H = [cos(2*pi*f0*t) sin(2*pi*f0*t)];
    size(C_w)
    
    theta_hat_gen_eq = al_est_bayes_mmse_linear_model(x,H,mu_theta,C_theta,C_w)
end

%% Sequential Estimator - Kay, Estimation Theory, example 12.3, p.399
% Sequential Estimator Initialization
theta_hat = zeros(2,1);
M = eye(2)* sigma_theta;
[theta_hat,theta_hat_v] = al_est_seq_lmmse(x,theta_hat,M,H,sigma_theta);

% I = eye(size(M,1));
% theta_hat_v = zeros(2,N);
% for k = 1:N
%     theta_hat_v(:,k) = theta_hat;
%     h = H(k,:)';
% %     [size(M) size(h)]
%     K = (M*h)/(sigma_theta.^2 + h'*M*h);
%     M = (I-K*h')*M;
%     theta_hat = theta_hat + K*(x(k)-h'*theta_hat);
% end

theta_hat
theta_hat_bayes_mmse
figure
plot(theta_hat_v'),grid on,hold on
return
