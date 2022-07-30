script
close all
clear global
clc

% rng(123)


%% generate GMM data
N = 1e5;
epsln = 0.05;
sigma1 = 1;
sigma2 = 30;
[x_gmm] = al_gen_gmm(N,epsln,sigma1,sigma2);

% figure
% plot(x_gmm,'b- .'),grid on,hold on


%% Estimate GMM epsilon using method of moments
epsln_hat = al_gmm_eps_est(x_gmm,sigma1,sigma2);
epsln
err = epsln_hat - epsln