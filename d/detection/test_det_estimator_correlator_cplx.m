%script
close all
clear
clc

% Test for Complex Estimator-Correlator Detector.
% Reference: 
% Kay, Fundamentals of Statistical Signal Processing, 
% Volume III Practical Algorithm Development,
% Algorithm 12.17 – Estimator-Correlator Detector (also Algorithm 10.4)

%% Generate random signal
Fs = 1;
N = 2^6;
sigma_n = 5.0;
var_wgn = sigma_n.^2

% Generate Correlated Random Signal
h = ones(8,1);
Nlags = length(h)*1-1;
R_h = conv(h,h);
R_h_1s = R_h(Nlags+1:end);
r = zeros(N,1);
r(1:length(R_h_1s)) = R_h_1s;
Cs = toeplitz(r);
% C = toeplitz(R_h(Nlags+1:end));
[s] = al_gen_corr_cwgn(Cs,1);

% Uncorrelated Noise Generator
w = randn(N,1)*sigma_n/sqrt(2) + 1j*randn(N,1)*sigma_n/sqrt(2);

% Add noise to signal
x = s + w;

% [norm(s) norm(w) norm(real(w)) norm(imag(w))]
% SNR_hat2 = 20*log10(norm(s)/norm(w))


%% Complex Estimator-Correlator Detector
% Pfa = 5e-2;
[Tx,s_hat] = al_det_estimator_correlator_cplx(x,Cs,var_wgn);

figure
plot(real(x),'r- .'),hold on,grid on
plot(real(s),'b- o'),hold on,grid on
plot(real(s_hat),'k-x'),hold on,grid on
legend({'x','s','s\_hat'})
title('Real part of Tx signal s and estimated signal s_hat')
figure
plot(imag(x),'r- .'),hold on,grid on
plot(imag(s),'b- o'),hold on,grid on
plot(imag(s_hat),'k-x'),hold on,grid on
legend({'x','s','s\_hat'})
title('Image part of Tx signal s and estimated signal s_hat')

% return

Nexp = 1e4;
Nfa = 0;% Number of False Alarm Events
Tx_v = zeros(Nexp,1);
Ts_v = zeros(Nexp,1);
Tn_v = zeros(Nexp,1);
thresh_v = zeros(Nexp,1);
for k = 1:Nexp
    w = randn(N,1)*sigma_n/sqrt(2) + 1j*randn(N,1)*sigma_n/sqrt(2);
    xn = w;
    xs = s + w;
    [Tn] = al_det_estimator_correlator_cplx(xn,Cs,var_wgn);
    [Ts] = al_det_estimator_correlator_cplx(xs,Cs,var_wgn);
    Ts_v(k) = Ts;
    Tn_v(k) = Tn;
%     thresh_v(k) = thresh;
%     if(Tx > thresh)
%         Nfa = Nfa + 1;
%     end
end
figure
plot(Tn_v,'r- .'),grid on,hold on
plot(Ts_v,'b- .'),grid on,hold on
legend({'H0','H1'})
title('Detector Statisics for H0 and H1')

[hist_Tn_v,x_Tn_v] = hist(Tn_v,100);
[hist_Ts_v,x_Ts_v] = hist(Ts_v,100);
figure
plot(x_Tn_v,hist_Tn_v,'r- .'),grid on,hold on
plot(x_Ts_v,hist_Ts_v,'b- .'),grid on,hold on
legend({'H0','H1'})
title('Histograms for H0 and H1')

%Compare given Pfa with Monte Carlo Test results
% Pfa_hat = Nfa/Nexp
% Pfa
% [mean(Tx_v) mean(thresh_v)]
% mean(Tx_v)/mean(thresh_v)

% p = (0:0.01:1);
% qi = al_q_inv_func(p);
% figure
% plot(p,qi,'b- .'),grid on,hold on

return
