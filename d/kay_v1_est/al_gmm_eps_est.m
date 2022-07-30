function epsln_hat = al_gmm_eps_est(x_gmm,sigma1,sigma2)
% Estimate epsilon parameter of the Gaussian Mixture Model
N = length(x_gmm);
x2 = x_gmm.^2;
sx2 = sum(x2);
epsln_hat = (sx2/N - sigma1^2)/(sigma2^2 - sigma1^2);
return