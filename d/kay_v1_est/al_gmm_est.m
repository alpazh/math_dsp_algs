function [epsln_hat,sigma1_sq_hat,sigma2_sq_hat] = al_gmm_est(x_gmm)
% Estimate epsilon parameter and sigma1, sigma2 of the Gaussian Mixture Model

mu2 = mean(x_gmm.^2);
mu4 = mean(x_gmm.^4);
mu6 = mean(x_gmm.^6);

u = (mu6-5*mu4*mu2)/(5*mu4-15*mu2*mu2);
v = mu2*u - mu4/3;
sigma1_sq_hat = (u + sqrt(u*u-4*v))/2;
sigma2_sq_hat = v/sigma1_sq_hat;
epsln_hat = 1 - (mu2 - sigma1_sq_hat)/(sigma2_sq_hat - sigma1_sq_hat);

return
