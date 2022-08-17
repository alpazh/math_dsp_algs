function [w] = al_gen_corr_cwgn(C,N)
% [w] = al_gen_corr_cwgn(C,N)
% Generate Complex Correlated (Colored) White Gaussian Noise
% provided Covariance Matrix C and number of correlated noise vectors N
% Input arguments:
% C - noise covariance matrix
% N - number of correlated noise vectors
% 
% Output arguments:
% w - generated noise
% 
% Author: Alexey Zherebtsov

% calculate number of dimensions M of vector to generate
M = size(C,1);
% generate N samples of M-dimensional random vectors
u = randn(M,N) + randn(M,N)*1j;
% make Choletsky factorization of the covariance matrix
F = chol(C);
F = F.';
% compute correlated random vectors
w = F*u;

end