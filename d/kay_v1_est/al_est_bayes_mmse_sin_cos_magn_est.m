function theta_hat = al_est_bayes_mmse_sin_cos_magn_est(x,sigma,sigma_theta,f0,Fs)
% Bayes MMSE estimator of the magnitude of sum sin+cos signals 
% of the same frequency in the AWGN.
% Kay, Estimation Theory, Chapter 11, p.347, Example 11.1:
N = length(x);
t = (0:1/Fs:(N-1)/Fs)';
H = [cos(2*pi*f0*t) sin(2*pi*f0*t)];
theta_hat = (1/sigma^2)/((1/sigma_theta^2)+(N/(2*sigma^2)))*(H'*x);
return
