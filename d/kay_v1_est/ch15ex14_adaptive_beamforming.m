script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 15, p.544, 
% Example 15.14 Adaptive Beamforming

M = 32;%number of sensors
F0 = 3e8;%carrier frequency
d = 0.125;%sensors spacing
c = 3e8;%propagation speed
beta = pi/3;%angle of arrival
phi = 0;%initial phase
t0 = 0;%initial time
A = 1.234;%signal amplitude

phi_prime = -2*pi*F0*t0+phi;
A_tilde = A*exp(1j*phi_prime);

fs = F0*(d/c)*cos(beta);%spatial frequency
e = exp(1j*2*pi*fs*(0:M-1)');

t = 0;
%noiseless signal
s_tilde_t = A_tilde*exp(1j*2*pi*F0*t)*e;
%noise

% define number of dimensions M of vector to generate
% M = 3
% define covariance matrix
% CM = randn(M,M)*diag(rand(M,1))/10
CM = (cos(2*pi*(1/128)*(1:M)+pi/16)/2);
C = CM'*CM;
[V,D] = eig(C);
C = abs(V*diag((M:-1:1)'/(M*2))*V');
% [V,D] = eig(C)
% return
% define signal length
N = 1;%3e3
w_tilde_t = al_gen_corr_cwgn(C,N);

% uncorrelated noise variant
% sigma_w = 1;
% C = eye(M)*sigma_w;
% w_tilde_t = (randn(M,1) + 1j*randn(M,1))*sigma_w;

%signal vector snapshot at the M sensors
x_tilde_t = s_tilde_t + w_tilde_t;

figure
imagesc(C)
title('Covariance Matrix')
figure
surfc(C)
title('Covariance Matrix')
figure
plot(real(s_tilde_t),'b- .'),grid on,hold on
plot(imag(s_tilde_t),'r- .'),grid on,hold off
figure
plot(real(x_tilde_t),'b- .'),grid on,hold on
plot(imag(x_tilde_t),'r- .'),grid on,hold off
title('Sensors output')

M = 32;%number of sensors
F0 = 3e8;%carrier frequency
d = 0.125;%sensors spacing
c = 3e8;%propagation speed
beta = pi/3;%angle of arrival
phi = 0;%initial phase
t0 = 0;%initial time
A = 1.234;%signal amplitude

%optimal beamformer
[y_tilde_t,a_opt] = al_mvdr_beamformer(x_tilde_t,C,F0,d,c,beta);

figure
plot(real(a_opt),'b- .'),grid on,hold on
plot(imag(a_opt),'r- .'),grid on,hold off
title('Beamformer weights')
figure
plot(real(y_tilde_t),'b- .'),grid on,hold on
plot(imag(y_tilde_t),'r- .'),grid on,hold off
title('Beamformer output')

return
