function A_hat_mmse = al_est_bayes_mmse_dc_level_in_wgn(x,sigma,sigma_A,mu_A)
% Ref: Kay, Estimation Theory, p.319, eq. 10.11
N = length(x);
alpha = sigma_A^2 / (sigma_A^2 + sigma^2/N);
A_hat_mmse = alpha*mean(x) + (1-alpha)*mu_A;
return
