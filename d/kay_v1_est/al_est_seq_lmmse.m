function [theta_hat,theta_hat_v] = al_est_seq_lmmse(x,theta_hat,M,H,sigma_theta)
% [theta_hat,theta_hat_v] = al_est_seq_lmmse(x,theta_hat,M,H,sigma_theta)  
% Estimates vector parameter theta_hat using sequential linear MMSE
% Estimator.
% Reference: Kay, Estimation Theory, Chapter 12.6 Sequential LMMSE
% Estimation.

%   Author: Alexey Zherebtsov

N = length(x);
I = eye(size(M,1));
theta_hat_v = zeros(2,N);
for k = 1:N
    theta_hat_v(:,k) = theta_hat;
    h = H(k,:)';
%     [size(M) size(h)]
    K = (M*h)/(sigma_theta.^2 + h'*M*h);
    M = (I-K*h')*M;
    theta_hat = theta_hat + K*(x(k)-h'*theta_hat);
end

return