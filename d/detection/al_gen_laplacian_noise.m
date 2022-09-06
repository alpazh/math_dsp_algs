function [w] = al_gen_laplacian_noise(var_lapl,N)
% w = al_gen_laplacian_noise(var_lapl,N)
% Laplacian Noise Generator
% It generates Laplacian Noise with given variance var_lapl.
%
% Input paramaters:
% var_lapl - Laplacian Noise Variance,
% N - number of samples to be generated.
%
% Output paramaters:
% w - generated noise.
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development,
% Exercise 4.13.

u = rand(N,1);
w = zeros(N,1);
w(u>0.5,1) = sqrt(var_lapl)*(1/sqrt(2))*log(1./(2*(1-u(u>0.5))));
w(u<=0.5,1) = sqrt(var_lapl)*(1/sqrt(2))*log(2*u(u<=0.5));

return
