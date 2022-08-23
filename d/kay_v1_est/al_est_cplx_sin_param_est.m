function [f0_hat,A_hat,phi_hat] = al_est_cplx_sin_param_est(x_n,Fs)
% [f0_hat,A_hat,phi_hat] = al_est_cplx_sin_param_est(x_n,f0)
% ML Sinusoidal Parameter Estimation in CWGN
% 
% Input arguments:
% x_n - input signal
% 
% Output arguments:
% f0_hat  - estimated frequency
% A_hat   - estimated magnitude
% phi_hat - estimated phase
% 
% Reference: Kay, Estimation Theory, Chapter 15, p.531, 
% example 15.10 phase of complex sinusoid in CWGN
% 
% Author: Alexey Zherebtsov

N = length(x_n);

% Estimate f0 by periodogram
AX2 = abs(fft(x_n)).^2;
[~,imax_AX2] = max(AX2);
f0_hat = (imax_AX2-1)*Fs/N;

% figure
% plot(AX2,'r- .'),grid on,hold on

% Estimate A
N = length(x_n);
t_v = (1/Fs)*(0:N-1)';
A_tilde_hat = (1/N)*(x_n.'*exp(-1j*2*pi*f0_hat*t_v));
A_hat = abs(A_tilde_hat);

% Estimate phi
phi_hat = atan2(imag(A_tilde_hat),real(A_tilde_hat));

end
