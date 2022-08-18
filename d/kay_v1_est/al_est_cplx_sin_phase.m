function [phi_hat] = al_est_cplx_sin_phase(x_n,f0)
% [phi_hat] = al_est_cplx_sin_phase(x_n,f0)
% ML Estimator of the phase of complex sinusoid in CWGN
% 
% Input arguments:
% x_n - input signal
% f0 - sinusoid frequency
% 
% Output arguments:
% phi_hat- estimated phase
% 
% Reference: Kay, Estimation Theory, Chapter 15, p.531, 
% example 15.10 phase of complex sinusoid in CWGN
% 
% Author: Alexey Zherebtsov

N = length(x_n);
t_v = (0:N-1)';
X_f0 = x_n.'*exp(-1j*2*pi*f0*t_v);
phi_hat = atan2(imag(X_f0),real(X_f0));

end
