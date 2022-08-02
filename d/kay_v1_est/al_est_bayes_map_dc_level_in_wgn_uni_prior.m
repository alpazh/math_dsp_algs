function A_hat_map = al_est_bayes_map_dc_level_in_wgn_uni_prior(x,A0)
% Estimate DC level in AWGN assuming uniform prior PDF.
% Kay, Estimation Theory, Chapter 11, p.355, Example 11.3

A_hat_map = mean(x);
if A_hat_map > A0
    A_hat_map = A0;
end
if A_hat_map < -A0
    A_hat_map = -A0;
end
return
