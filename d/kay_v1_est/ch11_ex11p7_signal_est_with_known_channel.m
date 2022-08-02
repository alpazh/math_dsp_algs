script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 11.7, p.367:
% MMSE Signal Estimation from the output of the known channel

% rng(123)
ns = 128;
A =  1.0;
sigma = 0.5;

% s = ((randn(N,1) > 0) - 0.5)*2*A;
% s0 = hamming(ns);
s0 = filter(ones(8,1)/8,1,randn(ns,1));
s0 = s0(1:ns,1);

% s0 = (1:ns)'/ns;
% s0 = randn(ns,1);

% h = (1:Nch/2)'*2/Nch;
% h = [h;flipud(h)];

% s = [zeros(ns,1);s0;zeros(ns,1)];
% s = [s0;zeros(ns,1)];
s = s0;
ns = length(s)

nh = 32
% h = zeros(ns,1);
h = (1:nh)'/nh;
% h(1:length(s0)) = flipud(s0);
% h = hamming(nh);
% nh = length(h)

N = ns + nh - 1

r = [h(1) zeros(1,ns-1)];
c = zeros(N,1);
c(1:nh,1) = h;
H  = toeplitz(c,r);

sh = conv(s,h);
x = sh + randn(length(sh),1)*sigma;
% x = x(1:ns,1);

figure
plot(s,'b- .'),grid on,hold on
figure
plot(h,'k- .'),grid on,hold on
figure
imagesc(H)
figure
plot(sh,'r- .'),grid on,hold on
% figure
plot(H*s,'m- o'),grid on,hold on

size_H = size(H)
size_Ht = size(H')
acs = xcorr(s);
acs = acs((length(acs)+1)/2:end);
Cs = toeplitz(acs);
figure
plot(acs,'r- .'),grid on
figure
imagesc(Cs)

s_hat = al_mmse_chan_est(Cs,H,sigma,x);
r_trunc = 30;
s_hat_tr = al_mmse_chan_est_trunc_svd(Cs,H,x,r_trunc);
% size_HxCsxHt = size(H*Cs*H')
% I = eye(N);
% size(x)
% [size(H*Cs*H')  size(I)]
% B = inv(H*Cs*H'+sigma^2*I);
% [size(Cs*H') size(B) size(x)]
% MTI = H*Cs*H'+sigma^2*I;
% cond_MTI = cond(MTI)
% s_hat = Cs*H'*inv(H*Cs*H'+sigma^2*I)*x;
% % Cs*H'*inv([ns ns]*Cs*H'+sigma^2*I)*[ns 1];
% size(s_hat)

figure
plot(x,'k- .'),grid on,hold on
figure
plot(s,'b- o'),grid on,hold on
plot(s_hat,'r- x'),grid on,hold on
figure
plot(s,'b- o'),grid on,hold on
plot(s_hat_tr,'m- d'),grid on,hold on

return
