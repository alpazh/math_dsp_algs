%script
close all
clear
clc

% Test for Replica-correlator (matched filter) detector.
% Reference: 
% Kay, Fundamentals of Statistical Signal Processing, 
% Volume III Practical Algorithm Development,
% Algorithm 12.14 – Replica-correlator (matched filter) (also Algorithm 10.1)

%% Generate signal
Fs = 1;
N = 2^6;
sigma_n = 3;
F1 = Fs/N;
t =(0:N-1)';
A = 1;
phi = pi*5/3;
s = A*exp(1j*(2*pi*F1*t+phi));
w = randn(N,1)*sigma_n/sqrt(2) + 1j*randn(N,1)*sigma_n/sqrt(2);
x = s + w;

% [norm(s) norm(w) norm(real(w)) norm(imag(w))]
% SNR_hat2 = 20*log10(norm(s)/norm(w))


%% Replica-Correlator Detector
var_wgn = sigma_n.^2
% w'*w
% return
Pfa = 5e-2;
[Tx,thresh] = al_det_replica_correlator_cplx(s,x,var_wgn,Pfa)

figure
plot(real(x),'b- .'),hold on,grid on
plot(imag(x),'r- .'),hold on,grid on
plot(real(s),'c-'),hold on,grid on
plot(imag(s),'m-'),hold on,grid on
title('Rx signal')

Nexp = 1e5;
n_sc_f = sigma_n/sqrt(2);%noise_scale_factor
Nfa = 0;% Number of False Alarm Events
Tx_v = zeros(Nexp,1);
thresh_v = zeros(Nexp,1);
for k = 1:Nexp
    w = (randn(N,1) + 1j*randn(N,1))*n_sc_f;
    x = w;
    [Tx,thresh] = al_det_replica_correlator_cplx(s,x,var_wgn,Pfa);
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
