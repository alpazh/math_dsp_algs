%script
close all
clear
clc

% Test for Limiter-Replica-Correlator (matched filter) detector.
% Reference: 
% Kay, Fundamentals of Statistical Signal Processing, 
% Volume III Practical Algorithm Development,
% Algorithm 12.15 – Limiter-Replica-Correlator (real-valued version is Algorithm 10.2)

%% Generate signal
Fs = 1;
N = 2^7;
var_lapl = 0.01;
% F1 = Fs/N;
t =(0:N-1)';
A = 1;
% phi = pi*5/3;
s = A*ones(N,1);
% Generate Laplacian Noise
% w = randn(N,1)*sigma_n;
[w] = al_gen_laplacian_noise(var_lapl,N);
% size(w)
% figure
% plot(w,'r- .'),grid on,hold on
% return
x = s + w;
% [norm(s) norm(w) norm(real(w)) norm(imag(w))]
% SNR_hat2 = 20*log10(norm(s)/norm(w))


%% Replica-Correlator Detector
Pfa = 5e-2;
[Tx,thresh] = al_det_limiter_replica_correlator(A,x,var_lapl,Pfa);

figure
plot(x,'b- .'),hold on,grid on
plot(s,'c-'),hold on,grid on
title('Rx signal')

% return

Nexp = 1e4;
Nfa = 0;% Number of False Alarm Events
Tx_v = zeros(Nexp,1);
thresh_v = zeros(Nexp,1);
for k = 1:Nexp
    [w] = al_gen_laplacian_noise(var_lapl,N);
    x = w;
    [Tx,thresh] = al_det_limiter_replica_correlator(A,x,var_lapl,Pfa);
    Tx_v(k) = Tx;
    thresh_v(k) = thresh;
    if(Tx > thresh)
        Nfa = Nfa + 1;
    end
end
figure
plot(Tx_v,'r- .'),grid on,hold on
plot(thresh_v.^1,'b- .'),grid on,hold on
Pfa_hat = Nfa/Nexp
Pfa
% [mean(Tx_v) mean(thresh_v)]
% mean(Tx_v)/mean(thresh_v)

% p = (0:0.01:1);
% qi = al_q_inv_func(p);
% figure
% plot(p,qi,'b- .'),grid on,hold on

return
