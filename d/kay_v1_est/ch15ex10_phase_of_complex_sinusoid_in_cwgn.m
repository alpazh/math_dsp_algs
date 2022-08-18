script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 15, p.531, 
% Example 15.10 Phase of complex sinusoid in CWGN

N = 2^8
% signal generation 
phi = 1*pi/3
A = 12.3
f0 = 1/(2^5)
w_n = randn(N,1)+randn(N,1)*1j;
t_v = (0:N-1)';
x_n = A*exp(1j*(2*pi*f0*t_v+phi)) + w_n;

figure
plot(real(x_n),'b- .'),grid on,hold on
plot(imag(x_n),'r- .'),grid on,hold off

% phase estimator
[phi_hat] = al_est_cplx_sin_phase(x_n,f0);

phi_hat/pi
return
