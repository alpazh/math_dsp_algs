script
close all
clear global
clc

% rng(123)


%% generate GMM data
N = 1e6;
epsln = 0.75;
sigma1 = 1;
sigma2 = 10;
[x_gmm] = al_gen_gmm(N,epsln,sigma1,sigma2);

% figure
% plot(x_gmm,'b- .'),grid on,hold on


%% Estimate GMM epsilon using method of moments
[epsln_hat,sigma1_sq_hat,sigma2_sq_hat] = al_gmm_est(x_gmm)
epsln_hat
epsln
err = epsln_hat - epsln