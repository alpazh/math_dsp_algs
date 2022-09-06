%script
close all
clear
clc

% Test for Generalized Matched Filter Detector.
% Reference: 
% Kay, Fundamentals of Statistical Signal Processing, 
% Volume III Practical Algorithm Development,
% Algorithm 12.16 – Generalized Matched Filter Detector (also Algorithm 10.3)

%% Generate signal
Fs = 1;
N = 2^6;
sigma_n = 0.5;
var_wgn = sigma_n.^2
F1 = Fs/N;
t =(0:N-1)';
A = 1;
phi = pi*5/3;
s = A*exp(1j*(2*pi*F1*t+phi));

% h = [0.9 0.8 0.6 0.3]';
h = ones(5,1);
Nlags = length(h)*1-1;

% Correlated Noise Generator Variant I: FIR
% v = randn(N,1)*sigma_n/sqrt(2) + 1j*randn(N,1)*sigma_n/sqrt(2);
% w = filter(h,1,v);
% r_ww = xcorr(w,Nlags);
% C_hat = toeplitz(r_ww(Nlags+1:end));
% figure
% plot(abs(r_ww),'m- s'),grid on,hold on
% figure
% imagesc(abs(C_hat))
% R_hat = real(C_hat)/(sigma_n.^2*N)

% Correlated Noise Generator Variant II: Choletsky
R_h = conv(h,h);
R_h_1s = R_h(Nlags+1:end);
r = zeros(N,1);
r(1:length(R_h_1s)) = R_h_1s;
C = toeplitz(r);
% C = toeplitz(R_h(Nlags+1:end));
[w] = al_gen_corr_cwgn(C,1);
% [size(s) size(w)]
% r_ww = xcorr(w,Nlags*2);
% figure
% plot(abs(r_ww),'m- s'),grid on,hold on
% title('Correlated Noise ACF')

x = s + w;

% [norm(s) norm(w) norm(real(w)) norm(imag(w))]
% SNR_hat2 = 20*log10(norm(s)/norm(w))


%% Generalized Matched Filter Detector

% w'*w
% return
Pfa = 5e-2;
[Tx,thresh] = al_det_generalized_matched_filter_cplx(s,x,Pfa,C);

figure
plot(real(x),'b- .'),hold on,grid on
plot(imag(x),'r- .'),hold on,grid on
plot(real(s),'c-'),hold on,grid on
plot(imag(s),'m-'),hold on,grid on
title('Rx signal')

% return

Nexp = 5e4;
Nfa = 0;% Number of False Alarm Events
Tx_v = zeros(Nexp,1);
thresh_v = zeros(Nexp,1);
for k = 1:Nexp
%     w = (randn(N,1) + 1j*randn(N,1))*n_sc_f;
    w = al_gen_corr_cwgn(C,1);
    x = w;
    [Tx,thresh] = al_det_generalized_matched_filter_cplx(s,x,Pfa,C);
    Tx_v(k) = Tx;
    thresh_v(k) = thresh;
    if(Tx > thresh)
        Nfa = Nfa + 1;
    end
end
figure
plot(Tx_v,'r- .'),grid on,hold on
plot(thresh_v,'b- .'),grid on,hold on

%Compare given Pfa with Monte Carlo Test results
Pfa_hat = Nfa/Nexp
Pfa
% [mean(Tx_v) mean(thresh_v)]
% mean(Tx_v)/mean(thresh_v)

% p = (0:0.01:1);
% qi = al_q_inv_func(p);
% figure
% plot(p,qi,'b- .'),grid on,hold on

return
