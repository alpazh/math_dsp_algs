%script
close all
clear
clc

% Test for GLRT Detector which detects complex exponent with known frequency and 
% unknown amplitude.
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development,
% Algorithm 12.20 – Unknown complex amplitude of complex exponent (also Algorithm 10.6)

%% Generate signal
Fs = 1;
N = 2^6;
sigma_n = 3;
var_wgn = sigma_n.^2;
F1 = 2*Fs/N;
t =(0:N-1)';
A = ceil(rand(1,1)*4+1)
phi = pi*5/3;
s = A*exp(1j*(2*pi*F1*t+phi));
w = randn(N,1)*sigma_n/sqrt(2) + 1j*randn(N,1)*sigma_n/sqrt(2);
x = s + w;

% [norm(s) norm(w) norm(real(w)) norm(imag(w))]
% SNR_hat2 = 20*log10(norm(s)/norm(w))

%% GLRT Detector
Pfa = 5e-2;
[Tx,thresh] = al_det_glrt_unk_ampl_unk_freq_cplx_exp(x,var_wgn,Pfa)

figure
plot(real(x),'r- .'),hold on,grid on
plot(real(s),'b- o'),hold on,grid on
legend({'x','s'})
title('Real part of Tx signal s and Rx signal x')
figure
plot(imag(x),'r- .'),hold on,grid on
plot(imag(s),'b- o'),hold on,grid on
legend({'x','s'})
title('Image part of Tx signal s and Rx signal x')

% return

Nexp = 1e5;
Nfa = 0;% Number of False Alarm Events
Nfa_k = 0;% Number of False Alarm Events
Tx_v = zeros(Nexp,1);
Ts_v = zeros(Nexp,1);
Tn_v = zeros(Nexp,1);
thresh_v = zeros(Nexp,1);
for k = 1:Nexp
    w = randn(N,1)*sigma_n/sqrt(2) + 1j*randn(N,1)*sigma_n/sqrt(2);    
    xn = w;
    xs = s + w;
    [Tn,thresh_H0] = al_det_glrt_unk_ampl_unk_freq_cplx_exp(xn,var_wgn,Pfa);
    [Ts,thresh_H1] = al_det_glrt_unk_ampl_unk_freq_cplx_exp(xs,var_wgn,Pfa);
    Ts_v(k) = Ts;
    Tn_v(k) = Tn;
    thresh_v(k) = thresh_H0;
    if(Tn > thresh_H0)
        Nfa = Nfa + 1;
    end
end
figure
plot(Tn_v,'r- .'),grid on,hold on
% plot(Ts_v,'b- .'),grid on,hold on
plot(thresh_v,'k- .'),grid on,hold on
legend({'H0','Threshold for H0'})
% legend({'H0','H1','Threshold for H0'})
title('Detector Statisics for H0 and H1')

[hist_Tn_v,x_Tn_v] = hist(Tn_v,100);
[hist_Ts_v,x_Ts_v] = hist(Ts_v,100);
figure
plot(x_Tn_v,hist_Tn_v,'r- .'),grid on,hold on
plot(x_Ts_v,hist_Ts_v,'b- .'),grid on,hold on
legend({'H0','H1'})
title('Histograms for H0 and H1')

%Compare given Pfa with Monte Carlo Test results
disp('Compare estimated Pfa_hat and given Pfa')
disp('to make sure that detectors works as expected:')
Pfa_hat = Nfa/Nexp
Pfa
% [mean(Tx_v) mean(thresh_v)]
% mean(Tx_v)/mean(thresh_v)

% p = (0:0.01:1);
% qi = al_q_inv_func(p);
% figure
% plot(p,qi,'b- .'),grid on,hold on

return
