script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 13, p.436, 
% Example 13.3 Time-Varying Channel Estimation based on Vector Kalman Filter

%% State model:
p = 2;% state model dimension
r = 2;% state noise dimension
M = 1;% obseravtion vector dimension
A0 = 1;
N = 10e1;
sigma_s = 1;
mu_u = 0;
mu_w = 0;
mu_s = A0;

sigma_u = (0.0001);
sigma_w = (0.01);
% sigma_n = (0.9).^((1:N)'+1);
% M_n_1_n_1 = 1;

A = diag([0.99 0.999]);%rand(p,p)/2
% A = diag([1 1 1]);
B = rand(p,r)/8;
H = zeros(M,p);%rand(M,p)/2
% H(1,1) = 1;H(2,2) = 1;
H(1,1) = 0.5;H(1,2) = -0.5;

v = (rem((0:N-1)',10)>4)*1.0;

figure
plot(v,'b-o'),grid on,hold on
title('Channel Input')

% return

uv = randn(r,N)*sigma_u;
Q = eye(r)*sigma_u;% Noise u Covariance 
wv = randn(M,N)*sigma_w;
C = eye(M)*sigma_w;% Observation Noise w Covariance 

sv = zeros(p,N);
xv = zeros(M,N);

% s(1) = randn(1,1)*sigma_s + mu_s;
% sv(:,1) = randn(p,1)*sigma_s + mu_s;
% h0 = [0.9 0.45]'
h0 = [1 0.9]'
sv(:,1) = h0; %initial channel coefficients
% x = zeros(N,1);

% sigma_n = (1/2)*ones(N,1);
% u = randn(N,1)*sigma_u + mu_u;
% w = randn(N,1).*sigma_n + mu_w;
% x(1) = s(1) + w(1);

H_arr = zeros(N,p);%time varying observation matrix is composed from 
% shifted input signal samples.

for n = 2:N
    s_n_1 = sv(:,n-1);
    u = uv(:,n);
    w = wv(:,n);
%     H = Hv(:,:,n)%Time-varying variant

    % State Propagation
    s_n = A*s_n_1 + B*u;

    H = v(n-1:n+p-2,1)';
    H_arr(n,:) = H;
    
    % Observation
    x_n = H*s_n + w;
        
    sv(:,n) = s_n;
    xv(:,n) = x_n;
end

figure
plot(sv'),grid on,hold on
title('p-th order Gauss-Markov model output')

figure
plot(xv'),grid on,hold on
title('Observation vector x')


%% Kalman Filter
M_n_1_n_1 = eye(p);
s_hat_n_1_n_1 = zeros(p,1);

[s_hat_n_n_v,M_n_n_v,M_n_n_1_v,MMSE_v,K_arr] =...
    al_vector_kalman_filter_chan_est(xv,A,B,H_arr,Q,C,M_n_1_n_1,s_hat_n_1_n_1);

% s_hat_n_n_v = zeros(p,N);
% M_n_n_v = zeros(p,p,N);
% M_n_n_1_v = zeros(p,p,N);
% MMSE_v = zeros(p,p,2*N);
% MMSE_v(:,:,2*1) = M_n_1_n_1;
% M_n_n_v(:,:,1) = M_n_1_n_1;
% 
% IM = eye(M);
% for k = 2:N
%     x_n = xv(:,k);
%     
%     % Prediction
%     s_hat_n_n_1 = A*s_hat_n_1_n_1;
%     
%     % Min Prediction MSE
%     M_n_n_1 = A*M_n_1_n_1*A' + B*Q*B';
%     M_n_n_1_v(:,:,k) = M_n_n_1;
%     MMSE_v(:,:,2*k-1) = M_n_n_1;
%     
%     % Kalman Gain
%     K_n = M_n_n_1*H'*inv(C + H*M_n_n_1*H');
%     
%     % Correction
%     s_hat_n_n = s_hat_n_n_1 + K_n * (x_n - H*s_hat_n_n_1);
%     s_hat_n_1_n_1 = s_hat_n_n;
%     
%     % Min MSE
% %     size(K_n*H)
% %     pause
%     M_n_n = (IM - K_n*H)*M_n_n_1;
%     M_n_1_n_1 = M_n_n;
%     
%     s_hat_n_n_v(:,k) = s_hat_n_n;
%     M_n_n_v(:,:,k) = M_n_n;    
%     MMSE_v(:,:,2*k) = M_n_n;
% end

figure
plot(sv','b- o'),grid on,hold on
plot(s_hat_n_n_v','r- x'),grid on,hold on
% legend({'s(n) to be estimated','estimated s(n)(Kalman Filter Output)'})
title('Kalman Filter Output')

% figure
% plot(log10(M_n_n_v),'r- ^'),grid on,hold on
% title('log10(Kalman Filter MSE)')
% 
% figure
% plot(M_n_n_1_v,'m- .'),grid on,hold on
% plot(M_n_n_v,'k- .'),grid on,hold on
% legend({'Min Prediction MSE','Min MSE'})
% title('Kalman Filter MSE')
% 
figure
mmse11 = MMSE_v(1,1,2:end);
% size(mmse11)
mmse11 = permute(mmse11,[3 2 1]);
% size(mmse11)
plot(mmse11','k- .'),grid on,hold on
title('Kalman Filter MSE[1,1]')

figure
plot(K_arr'),grid on,hold on
title('Kalman Filter Gains')

return
