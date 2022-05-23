close all
clear
clc
% fix the random seed to reproduce results
rng(465231)

%% Generate Gaussian random data for given covariance

% define number of dimensions M of vector to generate
M = 3
% define covariance matrix
CM = randn(M,M)*diag(rand(1,M));
CM = CM'*CM
% define number of data points
N = 3e4
% generate N samples of M-dimensional random vectors
x = randn(M,N);
% make Cholesky factorization of the covariance matrix
F = chol(CM)

% compute correlated random vectors
F = F.'
y = F*x;

figure
yt = y.';
plot3(yt(:,1),yt(:,2),yt(:,3),'b .'),grid on,hold on
title('Data with 3 independent experiments')

%% Estimate covariance matrix from generated data to compare with given cov
% estimate sample mean
mu_est = mean(y')
% remove mean from data
y_unb = y - mu_est';
CM_est = y_unb*y_unb'/(N-1)
CM

% The Covariance Matrix V is Positive Semidefinite
[V,D] = eig(CM_est);
% all eigenvalues are real and not negative
ev_cov = diag(D)

% The covariance matrix is positive definite 
% unless the experiments are dependent
yd = y;
yd(2,:) = yd(3,:) + randn(size(yd(3,:)))*0.01;

mu_est = mean(yd')
% remove mean from data
yd_unb = yd - mu_est';
CMd_est = yd_unb*yd_unb'/(N-1)

[Vd,Dd] = eig(CMd_est);
% all eigenvalues are real and not negative
% dependent experiments produces 0 eigenvalue
ev_cov_d = diag(Dd)

figure
yt = yd.';
plot3(yt(:,1),yt(:,2),yt(:,3),'b .'),grid on,hold on
title('Data when 2 out of 3 experiments are dependent')

%% Diagonalizing the covariance matrix V means finding M independent
% experiments as combinations of the original M experiments.
[V,D] = eig(CM_est);
CM_est_diag = V'*CM_est*V
% Decorrelate data
[size(V') size(y)]
y_dec = V'*y;
mu_est = mean(y_dec')
% remove mean from data
y_dec_unb = y_dec - mu_est';
CM_dec_est = y_dec_unb*y_dec_unb'/(N-1)

figure
yt = y_dec.';
plot3(yt(:,1),yt(:,2),yt(:,3),'b .'),grid on,hold on
title('Decorrelated Data')

%% The Mean and Variance of z = x + y
% Mean of sum = Sum of means
yb = y; 
yb(1,:) = yb(1,:)+3;
yb(2,:) = yb(2,:)+5;

m1 = mean(yb(1,:))
m2 = mean(yb(2,:))
m_1plus2 = mean(yb(1,:)+yb(2,:))

% The variance of z = x + y is sz^2 = sx^2 + sy^2 + 2s_xy.
sz2 = std(y(1,:)+y(2,:)).^2
sz2_from_cov_mat = CM_est(1,1)+CM_est(2,2)+2*CM_est(1,2)

%% The Covariance Matrix for Z = AY
% The covariance matrix of Z = A*Y is Vz = A*Vy*A'
A = rand(M,M)
z = A*y;
CM_z = A*CM*A'
% estimate z covariance 
% to demonstrate the correctness of CM_z
mu_est = mean(z');
z_unb = z - mu_est';
CM_z_est = z_unb*z_unb'/(N-1)

%% The Correlation p
% p_xy = s_xy/(s_x*s_y)
p_y1_y2 = CM(1,2)/sqrt(CM(1,1)*CM(2,2))
% Correlation p is cosine of angle between data vectors
p_y1_y2_est = (y(1,:)*y(2,:)')/sqrt((y(1,:)*y(1,:)')*(y(2,:)*y(2,:)'))

return
