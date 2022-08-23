script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 15, p.531, 
% Example 15.13 Sinusoidal Parameter Estimation

%% Signal generation 
N = 2^8
Fs = 1e6

% Sinusoidal Signal Parameters
phi = 2*pi/3
A = 3.21*exp(1j*phi)
f0 = Fs*1/(2^5)

sigma_w = 0.25;

w_n = (randn(N,1) + randn(N,1)*1j)*sigma_w;
t_v = (1/Fs)*(0:N-1)';
x_n = A*exp(1j*(2*pi*f0*t_v)) + w_n;

figure
plot(real(x_n),'b- .'),grid on,hold on
plot(imag(x_n),'r- .'),grid on,hold off

%% ML Sinusoidal Signal Parameters Estimator
[f0_hat,A_hat,phi_hat] = al_est_cplx_sin_param_est(x_n,Fs);
disp('Estimated Parameters')
disp('f0_hat and f0:')
[f0_hat f0]
disp('A_hat and abs(A):')
[A_hat abs(A)]
disp('phi_hat and phi:')
[phi_hat phi]

return
