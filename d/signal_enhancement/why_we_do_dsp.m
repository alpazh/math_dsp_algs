%Little demo of the DSP strength
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
L = 256
L_1 = L - 1

% the desired signal is a harmonic random process
% clean signal parameters and generation
Fs = 48e3;
Np = 100*Fs;
Np_to_plot = 100;
Np_to_sound = 3*Fs;

Ts = 1/Fs;
F0 = 5e3%(2/5)*Fs/2
A = 0.1
phi_x = rand(1,1)*2*pi
t = (Ts:Ts:Np*Ts)';

x = A*cos(2*pi*F0*t + phi_x);
% noise parameters and generation
Mu_n = 0;
Sigma_n = 0.2;
figure
plot(x(L+1:L+Np_to_plot),'k- .'),grid on, hold off
title('this is useful signal...')
disp('this is useful signal...')
pause(1)
sound(x(L+1:L+Np_to_sound),Fs)
disp('press any key to continue...')
pause


v = randn(Np,1)*Sigma_n + Mu_n;
figure
plot(v(L+1:L+Np_to_plot),'r- .'),grid on, hold off
title('and this is noise which will hide our signal...')
disp('and this is noise which will hide our signal...')
pause(1)
sound(v(L+1:L+Np_to_sound),Fs)
disp('press any key to continue...')
pause

% noise plus signal mixture
y = x + v;
figure
plot(y(L+1:L+Np_to_plot),'b- .'),grid on, hold on
title('noise hid our signal, can you still hear it or not?')
disp('noise hid our signal, can you still hear it or not?')
pause(1)
sound(y(L+1:L+Np_to_sound),Fs)
disp('press any key to continue...')
pause


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

%figure
%imagesc(Ry)
figure
imagesc(Rx_p_v)
title('This matrix will help us to find the signal')

%Ry./Rx_p_v

% Wiener filter
i_i = zeros(L,1);
i_i(1) = 1;
h_w = inv(Ry)*Rx*i_i;
x_est = filter(h_w,1,y);

figure
plot(x_est(L+1:L+Np_to_plot),'b- .'),grid on, hold on
title('To save our signal from the noise, we did some DSP and filtered noise out, can you hear the signal again?')
disp('To save our signal from the noise, we did some DSP and filtered noise out, can you hear the signal again?')
pause(1)
sound(x_est(L+1:L+Np_to_sound),Fs)
disp('press any key to continue...')
pause

%figure
%plot(x,'b- s'),grid on,hold on
%plot(x_est,'r- x'),grid on,hold off
SNRout_est = 20*log10(rms(x)/rms(x_est-x))

% Wiener filter another formulation
ccf_yx = xcorr(y,x,L_1);
p_xy = ccf_yx(L_1+1:end,1);
h_w2 = inv(Ry)*p_xy;
size(h_w)
size(h_w2)
% but same result
figure
plot(h_w,'b- s'),grid on,hold on
plot(h_w2,'k- x'),grid on,hold off
title('Filter coefficients. Do they look similar to something you saw before?')
disp('That''s all folks.')
