close all
clear
clc
% pkg load signal% comment out in MatLab

%% Upward shift cyclic permutation
N = 4
x = (1:N)'
I = eye(N)
% Upward cyclic shift matrix P
P = [I(:,N) I(:,1:N-1)]
x1 = P*x
x2 = P*P*x
x3 = P*P*P*x
x4 = P*P*P*P*x
% Upward cyclic shift matrix P Properties 
% P^N = I
P4 = P*P*P*P

%% Circulant Matrix
c = (1:N)'*1
d = (1:N)'*10
c = randn(N,1)
d = randn(N,1)
% Circulant Matrix C is combination of P P^2 ... P^(N-1)
C = c(1)*I+c(2)*P+c(3)*P*P+c(4)*P*P*P
D = d(1)*I+d(2)*P+d(3)*P*P+d(4)*P*P*P
% Circulant Matrix C Properties
% It has constant diagonals
figure
imagesc(C)
figure
imagesc(D)
% Circulant Matrices Product is Circulant Matrix
G = C*D
figure
imagesc(G)
% Circulant Matrices Product is cyclic convolution of vectors c and d
res_lin_conv = conv(c,d)
res_cycl_conv = cconv(c,d,N)

%% Eigenvalues and Eigenvectors of P
% The eigenvalues of P are the N-th roots of 1:
[V,L] = eig(P);
ev = diag(L)
figure
plot(real(ev),imag(ev),'b o')
% eigenvector matrix for P is the Fourier matrix
w = exp(2*pi*1i/N)
n = (0:N-1)';
nn = n*n'
% Columns of the Fourier matrix F are the same as 
% columns of eigenvectors matrix V (order could be different)
% Orthogonal matrices have orthogonal eigenvectors
F = (w.^nn)./sqrt(N)
Fi = inv(F)
% F*Fi
% return
V

%% Eigenvalues and Eigenvectors of C
[VC,LC] = eig(C);
[VD,LD] = eig(D);
% C Eigenvectors are the same as P eigenvectors.
VC
% The N eigenvalues of C are the components of 
% Fc = inverse Fourier transform of c
evc = diag(LC)

%% The Convolution Rule
% each eigenvalue of C*D is just 
% eigenvalue of C times eigenvalue of D
diag((Fi*C*F)*(Fi*D*F)).'
diag((Fi*C*D*F)).'

%% Cross-correlation and Autocorrelation
% Cross correlation could be computed by convolution.
CC_of_cd = xcorr(c,d)
CC_of_cd_conv = conv(c,flipud(d))
figure
plot(CC_of_cd,'b- o'),grid on,hold on
plot(CC_of_cd_conv,'r- x'),grid on,hold on
legend({'Cross-correlation','Cross-correlation by convolution'})
return