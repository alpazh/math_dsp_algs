function [y_tilde_t,a_opt] = al_mvdr_beamformer(x_tilde_t,C,F0,d,c,beta)
% [y_tilde_t,a_opt] = al_mvdr_beamformer(x_tilde_t,C,F0,d,c,beta)
% minimum variance distortionless response beamformer
% minimize a'*C*a subject to constraint e'*a = 1
% 
% Input arguments:
% x_tilde_t - signal vector snapshot at the sensors
% C - noise covariance matrix
% F0 - carrier frequency
% d - sensors spacing
% c - propagation speed
% beta - angle of arrival
% 
% Output arguments:
% y_tilde_t - beamformer output
% a_opt - beamformer weights
% 
% Reference: Kay, Estimation Theory, Chapter 15, 
% example 15.14 minimum variance distortionless response beamformer 
% 
% Author: Alexey Zherebtsov

M = length(x_tilde_t);% number of sensors

fs = F0*(d/c)*cos(beta);%spatial frequency
e = exp(1j*2*pi*fs*(0:M-1)');

% Ci = inv(C);
% Ci_e = inv(C)*e
Ci_e = C\e;
%beamformer optimal coefficients
a_opt = (Ci_e)/(e'*Ci_e);
%beamformer output
y_tilde_t = a_opt'*x_tilde_t;

end