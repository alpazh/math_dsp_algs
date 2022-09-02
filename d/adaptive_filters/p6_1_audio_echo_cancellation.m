script
close all
clear all
clc
% Adaptive Filters Ali H Sayed
% Project VI.1 Acoustic echo cancellation 

% load room impulse response
load('room_ir.mat','-ascii');
h = room_ir;
figure
plot(h,'b- '),grid on,hold on
title('Room Impulse Response')

% plot room frequency response
Fs = 8e3;
N = 1024;
[H,f] = freqz(h,1,N,Fs);

figure
plot(f,20*log10(abs(H)),'b- '),grid on,hold on
xlabel('frequency, Hz')
ylabel('20*log10(abs(H)), dB')
title('Frequency Response')

% load composite source signal
load('composite_source_signal.mat','-ascii');
% s = [composite_source_signal(1:2000,1); composite_source_signal(2801:4800,1) ];
s = composite_source_signal;

% repeat loudspeaker signal periodically
n_blocks = 24;
s = repmat(s,n_blocks,1);

s_echo = filter(h,1,s);

figure
plot(s,'k- '),grid on,hold on
title('Loudspeaker signal')

figure
plot(s_echo,'k- '),grid on,hold on
title('Echo signal')

size(s)
% sound(s)

%% Epsilon-NLMS AEC
M = 1024;% number of taps
w = zeros(M,1);% weights vector
mu = 1/(2.^0);% step size
d = s_echo;% desired signal is echo signal

[y,e,w] = al_eps_nlms(s,d,mu,w);

figure
plot(e,'b-'),grid on
title('suppressed echo signal')

figure
plot(w,'k-'),grid on
title('filter weights')

% sound(e)
return
