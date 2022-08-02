function [A_hat_map,sigma_2_hat] = al_est_bayes_map_dc_level_in_wgn_unkn_var(x,gamma,alpha,mu_A)
% Estimate DC level in AWGN with unknown variance assuming gamma prior PDF.
% Kay, Estimation Theory, Chapter 11, p.355, Example 11.4

N = length(x);
A_hat_map = (N*mean(x) + mu_A/alpha)/(N + 1/alpha);
sigma_2_hat = (N/(N+5))*(sum(x.^2)/N - A_hat_map^2) ...
    + (1/((N+5)*alpha))*(mu_A^2-A_hat_map^2)...
    + (2*gamma/(N+5));

return
