function theta_hat_mmse = al_est_bayes_mmse_linear_model(x,H,mu_theta,C_theta,C_w)
% Bayes LMMSE estimator for general linear model H*x+w=y
% Ref: Kay, Estimation Theory, ch.10, p.326, eq. 10.28

theta_hat_mmse =...
    mu_theta + C_theta*H'*(inv(H*C_theta*H'+C_w))*(x-H*mu_theta);

return
