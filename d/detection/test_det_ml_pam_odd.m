%script
close all
clear
clc

% PAM-M (for odd value of M) ML Detector
% This is Min prob of error (MPOE) detector
% which is reduced to MAP detector
% which is reduced to ML Detector
%
% Reference: Kay, Detection Theory, example 3.6

Nb = 1e5;
OSF = 5
N = Nb*OSF;
sigma_n = 0.5;
M = 3;


%% generate x (PAM-3 plus AWGN)
n = randn(N,1)*sigma_n;
s_symb = ceil(rand(Nb,1)*M);
s = s_symb - (M+1)/2;
s = repmat(s,1,OSF).';
s = s(:);
x = s + n;

SNR_hat2 = 20*log10(norm(s)/norm(n))

% calc s_mean for demo purposes
s_symb_resh = reshape(s,OSF,Nb);
s_mean = mean(s_symb_resh);

%% PAM-3 ML Detector 
% Min prob of error reduced to MAP reduced to ML Detector
% Assuming ideal synchronization
[s_symb_hat,x_mean] = al_det_ml_pam_odd(x,OSF,M);

Ns_to_plot = 1e2;
N_to_plot = Ns_to_plot*OSF;

figure
plot(s(1:N_to_plot),'b-o'),hold on,grid on
plot(x(1:N_to_plot),'r-x'),hold on,grid on
title('Tx and Rx signals')

figure
plot(s_mean(1:Ns_to_plot),'b-o'),hold on,grid on
plot(x_mean(1:Ns_to_plot),'r-x'),hold on,grid on
title('Estimated mean_x (per symbol)')

figure
plot(s_symb(1:Ns_to_plot),'b-s'),hold on,grid on
plot(s_symb_hat(1:Ns_to_plot),'r-x'),hold on,grid on
title('Detected symbol')

% size(s_symb)
% size(s_symb_hat)
% (s_symb~=s_symb_hat)
% nnz(s_symb~=s_symb_hat)
SER = nnz(s_symb~=s_symb_hat)/length(s_symb_hat)
return
