%Single-channel Signal Enhancement in the Time Domain
close all
clear
clc

if(exist('OCTAVE_VERSION', 'builtin'))
    pkg load communications
    pkg load signal
    pkg load statistics
end

% Signal Model and Problem Formulation
% filter length
L = 32
L_1 = L - 1

% plotting params
Nb = 2*L+1;
Len_to_plot = 60;

% the desired signal is a harmonic random process
% clean signal parameters and generation
Np = 5*48e3;
Fs = 48e3;
Ts = 1/Fs;
F0 = 4.8e3%(2/5)*Fs/2
A = 0.5
phi_x = rand(1,1)*2*pi
t = (Ts:Ts:Np*Ts)';

x = A*cos(2*pi*F0*t + phi_x);
% noise parameters and generation
Mu_n = 0;
Sigma_n = 1.25;

figure
plot(x(Nb+1:Nb+Len_to_plot),'k- .'),grid on, hold off
title('clean signal')

v = randn(Np,1)*Sigma_n + Mu_n;
figure
plot(v(Nb+1:Nb+Len_to_plot),'r- .'),grid on, hold off
title('noise')

% noise plus signal mixture
y = x + v;
figure
plot(y(Nb+1:Nb+Len_to_plot),'b- .'),grid on, hold on
title('noised signal')

% the input SNR
SNRin = 10*log10((A^2/2)/Sigma_n^2)
SNRin_est = 20*log10(rms(x)/rms(v))

% filter parameters estimation

% clean signal correlation matrix estimation
acf_x = xcorr(x,L_1);
acf_x_one_sided = acf_x(L_1+1:end,1);
Rx = toeplitz(acf_x_one_sided);
% noise correlation matrix estimation
acf_v = xcorr(v,L_1);
acf_v_one_sided = acf_v(L_1+1:end,1);
Rv = toeplitz(acf_v_one_sided);
% input signal correlation matrix estimation
acf_y = xcorr(y,L_1);
acf_y_one_sided = acf_y(L_1+1:end,1);

Ry = toeplitz(acf_y_one_sided);
Rx_p_v = Rx + Rv;

figure
imagesc(Ry)
title('Ry')

% pictures are about the same
% figure
% imagesc(Rx_p_v)
% title('Rx_p_v')

%Ry./Rx_p_v

% Wiener filter
i_i = zeros(L,1);
i_i(1) = 1;
h_w = inv(Ry)*Rx*i_i;
x_est = filter(h_w,1,y);

figure
plot(x_est(Nb+1:Nb+Len_to_plot),'b- .'),grid on, hold on
title('filtered signal')

figure
plot(x(Nb+1:Nb+Len_to_plot),'b- s'),grid on,hold on
plot(x_est(Nb+1:Nb+Len_to_plot),'r- x'),grid on,hold off
title('clean and filtered signal')

SNRout_est = 20*log10(rms(x)/rms(x_est-x))
% Y*H = X
% Y*Y*H = Y*X
% YY*H = YX
% H = inv(YY)*YX

% Wiener filter another formulation
ccf_yx = xcorr(y,x,L_1);
p_xy = ccf_yx(L_1+1:end,1);
%h_w2 = inv(Ry)*p_xy;
h_w2 = Ry\p_xy;
size(h_w)
size(h_w2)
% but same result
figure
plot(h_w,'b- s'),grid on,hold on
plot(h_w2,'k- x'),grid on,hold off
title('Wiener filter coefficients computed by two different ways')

% what about condition number of the matrix to be inverted?
Cond_num = cond(Ry)

% trade-off filter
% Lagrange Multiplier Value
% Mu = 0 is the identity filter
% Mu = 1 is the Wiener filter
% Mu > 1 is the filter with low  residual noise and high desired signal distortion
Mu = 4;
% Mu < 1 is the filter with high residual noise and low  desired signal distortion
% Mu = 1/100;
h_trade_off = inv(Rx + Mu*Rv)*Rx*i_i;
x_est_trade_off = filter(h_trade_off,1,y);

% MVDR filter
[V, Lambda] = eig(Rx);

lambda = diag(Lambda);

eig_val_max = max(lambda)
eig_val_thresh = eig_val_max*0.5

P = nnz(lambda > eig_val_thresh)
% Eigenvectors Matrix for nonzero eigenvalues
Qx_1 = V(:,lambda > eig_val_thresh);
size(Qx_1)

figure
plot(lambda),grid on
title('Rx eigenvalues')

figure
plot(Qx_1,'. -'),grid on
title('Rx eigenvectors corresponding to the signal subspace')

Rv_i = inv(Rv);
h_MVDR = Rv_i*Qx_1*inv(Qx_1'*Rv_i*Qx_1)*Qx_1'*i_i;
% h_MVDR can also be calcultated by another equation
Ry_i = inv(Ry);
h_MVDR2 = Ry_i*Qx_1*inv(Qx_1'*Ry_i*Qx_1)*Qx_1'*i_i;

% but results are the same
figure
plot(h_MVDR,'r- s'),grid on,hold on
plot(h_MVDR2,'b- x'),hold off
legend({'MVDR','MVDR2'})
title('MVDR filter coefficients computed by two different ways')

size(h_MVDR)
x_est_MVDR = filter(h_MVDR,1,y);

% TODO add maximum SNR filter
Rv_i_Rx = Rv_i*Rx;
[V2, Lambda2] = eig(Rv_i_Rx);
lambda2 = diag(Lambda2);
[Max_lambda2,Indx_lambda2_max] = max(lambda2)
t1 = V2(:,Indx_lambda2_max);
Zeta = 1/2; % any real value greater than zero
h_max_snr = Zeta*t1;
x_est_max_snr = filter(h_max_snr,1,y);

figure
plot(h_w,'r- s'),grid on,hold on
plot(h_trade_off,'g- o'),hold on
plot(h_MVDR,'b- d'),hold on
plot(h_max_snr,'k- x'),hold off
legend({'Wiener','Trade-off','MVDR','Max SNR'})
title('Filters coefficients comparison')

figure
plot(x_est(Nb+1:Nb+Len_to_plot),'r- s'),grid on,hold on
plot(x_est_trade_off(Nb+1:Nb+Len_to_plot),'g- o'),hold on
plot(x_est_MVDR(Nb+1:Nb+Len_to_plot),'b- d'),hold on
plot(x_est_max_snr(Nb+1:Nb+Len_to_plot),'k- x'),hold off
legend({'Wiener','Trade-off','MVDR','Max SNR'})
title('Filters output signals comparison')

