script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 13, p.420, 

%% 1st order Gauss-Markov model:
A0 = 5;
N = 1e2;
sigma_u = 0.01;
sigma_s = 0.1;
mu_u = 0;
mu_s = A0;
a = 0.999;

s = zeros(N,1);
s(1) = randn(1,1)*sigma_s + mu_s;
u = randn(N,1)*sigma_u + mu_u;
for k = 2:N
    s(k) = a*s(k-1) + u(k);
end

figure
plot(s,'b- .'),grid on,hold on
title('1st order Gauss-Markov model output')

%% p-th order State Model:
p = 5;
r = 3;
N = 2e1;
sigma_u = 0.01;
[~,U_s] = lu(rand(p,p))
C_s = U_s.'*U_s
UC_s = chol(C_s)
LS_s = UC_s.'

mu_u = 0;
mu_s = ones(p,1)*A0;

A = rand(p,p)/3
% [V,D] = eig(rand(p,p))
% A = real(V)*diag(rand(p,1))*real(V')

B = rand(p,r)/3
% [UB,SB,VB] = svd(rand(p,r))
% B = abs(UB*rand(p,r)*VB')

s_v = zeros(p,N);
s_1 = LS_s*randn(p,1) + mu_s;
u_v = randn(r,N)*sigma_u + mu_u;
s_v(:,1) = s_1;
for k = 2:N
    u = u_v(:,k);
    s = A*s_1 + B*u;
    s_1 = s;
    s_v(:,k) = s;
end

figure
plot(s_v'),grid on,title('p-th order Gauss-Markov model output')

return
