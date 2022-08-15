function [s_hat_n_n_v,M_n_n_v,M_n_n_1_v,MMSE_v] = ...
    al_vector_kalman_filter(xv,A,B,H,Q,C,M_n_1_n_1,s_hat_n_1_n_1)
% [s_hat_n_n_v,M_n_n_v,M_n_n_1_v,MMSE_v] = al_vector_kalman_filter(x,sigma_n,sigma_u,a,M_n_1_n_1,s_hat_n_1_n_1)
% Vector Kalman Filter
% 
% Input arguments:
% xv - observed vector signal 
% A - State Model Matrix A (State Propagation Matrix)
% B - State Model Matrix B (Driving Noise Matrix)
% H - Observation Matrix H
% Q - Driving Noise Covariance
% C - Observation Noise Covariance
% M_n_1_n_1 - initial MSE
% s_hat_n_1_n_1 - initial state vector estimate
% 
% Output arguments:
% s_hat_n_n_v - Estimated signal (Kalman Filter Output)
% M_n_n_v - Kalman Filter MSE
% M_n_n_1_v - Prediction MSE
% MMSE_v - vector of interleaved Kalman Filter MSE and Prediction MSE
% 
% Reference: Kay, Estimation Theory, Chapter 13.6 Vector Kalman Filter
% 
% Author: Alexey Zherebtsov

N = size(xv,2);
p = size(A,1);
M = size(H,1);
s_hat_n_n_v = zeros(p,N);
M_n_n_v = zeros(p,p,N);
M_n_n_1_v = zeros(p,p,N);
MMSE_v = zeros(p,p,2*N);
MMSE_v(:,:,2*1) = M_n_1_n_1;
M_n_n_v(:,:,1) = M_n_1_n_1;

IM = eye(M);
for k = 2:N
    x_n = xv(:,k);
    
    % Prediction
    s_hat_n_n_1 = A*s_hat_n_1_n_1;
    
    % Min Prediction MSE
    M_n_n_1 = A*M_n_1_n_1*A' + B*Q*B';
    M_n_n_1_v(:,:,k) = M_n_n_1;
    MMSE_v(:,:,2*k-1) = M_n_n_1;
    
    % Kalman Gain
%     K_n = (M_n_n_1*H')*inv(C + H*M_n_n_1*H');
%     Replace b*inv(A) with b/A
    K_n = (M_n_n_1*H')/(C + H*M_n_n_1*H');
    
    % Correction
    s_hat_n_n = s_hat_n_n_1 + K_n * (x_n - H*s_hat_n_n_1);
    s_hat_n_1_n_1 = s_hat_n_n;
    
    % Min MSE
    M_n_n = (IM - K_n*H)*M_n_n_1;
    M_n_1_n_1 = M_n_n;
    
    s_hat_n_n_v(:,k) = s_hat_n_n;
    M_n_n_v(:,:,k) = M_n_n;
    MMSE_v(:,:,2*k) = M_n_n;
end

end
