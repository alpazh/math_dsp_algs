script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 13, p.436, 
% Example 13.2 Scalar State Scalar Observation Kalman Filter

%% 1st order Gauss-Markov model:
A0 = 0;
N = 4e1;
sigma_s = 1;
mu_u = 0;
mu_w = 0;
mu_s = A0;

% Parameters for ex13.2:
% a = 1/2;
% sigma_u = sqrt(2);
% sigma_n = (1/2).^(1:N)';
% M_n_1_n_1 = 1;

% Parameters for Figure 13.7:
a = 0.99;
sigma_u = sqrt(0.1);
sigma_n = (0.9).^((1:N)'+1);
M_n_1_n_1 = 1;

s = zeros(N,1);
s(1) = randn(1,1)*sigma_s + mu_s;
x = zeros(N,1);

% sigma_n = (1/2)*ones(N,1);
u = randn(N,1)*sigma_u + mu_u;
w = randn(N,1).*sigma_n + mu_w;
x(1) = s(1) + w(1);

for k = 2:N
    % state propagation
    s(k) = a*s(k-1) + u(k);
    % observation
    x(k) = s(k) + w(k);  
end

figure
plot(s,'b- .'),grid on,hold on
title('1st order Gauss-Markov model output')

%% Kalman Filter
% M_n_1_n_1 = 1;
s_hat_n_1_n_1 = 0;

[s_hat_n_n_v,M_n_n_v,M_n_n_1_v,MMSE_v] = al_sc_st_sc_observ_kalman_f(x,sigma_n,sigma_u,a,M_n_1_n_1,s_hat_n_1_n_1);

% s_hat_n_n_v = zeros(N,1);
% M_n_n_v = zeros(N,1);
% M_n_n_1_v = zeros(N,1);
% MMSE_v = zeros(2*N,1);
% MMSE_v(2*1) = M_n_1_n_1;
% M_n_n_v(1) = M_n_1_n_1;
% for k = 2:N
%     % Prediction
%     s_hat_n_n_1 = a*s_hat_n_1_n_1;
%     % Min Prediction MSE
%     M_n_n_1 = a^2*M_n_1_n_1 + sigma_u^2;
%     M_n_n_1_v(k) = M_n_n_1;
%     MMSE_v(2*k-1) = M_n_n_1;
%     % Kalman Gain
%     K_n = M_n_n_1/(sigma_n(k)^2 + M_n_n_1);
%     % Correction
%     s_hat_n_n = s_hat_n_n_1 + K_n * (x(k) - s_hat_n_n_1);
%     s_hat_n_1_n_1 = s_hat_n_n;
%     % Min MSE
%     M_n_n = (1 - K_n)*M_n_n_1;
%     M_n_1_n_1 = M_n_n;
%     s_hat_n_n_v(k) = s_hat_n_n;
%     M_n_n_v(k) = M_n_n;    
%     MMSE_v(2*k) = M_n_n;
% end

figure
plot(s,'b- o'),grid on,hold on
plot(s_hat_n_n_v,'r- x'),grid on,hold on
legend({'s(n) to be estimated','estimated s(n)(Kalman Filter Output)'})
title('Kalman Filter Output')

figure
plot(log10(M_n_n_v),'r- ^'),grid on,hold on
title('log10(Kalman Filter MSE)')

figure
plot(M_n_n_1_v,'m- .'),grid on,hold on
plot(M_n_n_v,'k- .'),grid on,hold on
legend({'Min Prediction MSE','Min MSE'})
title('Kalman Filter MSE')

figure
plot(MMSE_v(2:end,1),'k- .'),grid on,hold on
title('Kalman Filter MSE')


return

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
plot(s_v'),grid on

return
