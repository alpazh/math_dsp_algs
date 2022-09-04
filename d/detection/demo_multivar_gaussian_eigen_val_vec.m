%script
close all
clear
clc

% Demo of Gaussian Multivariable Random Values 
% Covariance Matrix Eigenvalues and Eigenvectors

N = 1e6';
sigma_n = 1;
%% generate WSS process realization
n = randn(N,1)*sigma_n;
Nma = 16;
% h = ones(Nma,1)/Nma;
h = [1.0:-0.1:0.1]';
nc = filter(h,1,n);


%% Estimate correlation
Mlags = 2*Nma-1;
% Mlags = N;
r_nc = xcorr(nc,Mlags);
r_nc_one_sided = r_nc(Mlags+1:end,1);
if length(nc)<1e5
    figure
    plot(nc,'b- .'),grid on,hold on
end

figure
plot(r_nc,'r- .'),grid on,hold on
title('ACF')

figure
plot(r_nc_one_sided,'r- .'),grid on,hold on
title('ACF one-sided')

R = toeplitz(r_nc_one_sided.');

figure
imagesc(R)
% surfc(R)
title('Correlation matrix R')

% PSD
Pxx = fft(r_nc);

figure
plot(10*log10(abs(fftshift(Pxx))),'r- .'),grid on,hold on
title('log scale PSD')

% R*V=V*D
[Vc,D] = eig(R);

lambda_v = flipud(diag(D));
Vd = fliplr(Vc);

% figure
% plot(lambda_v,'g- .'),grid on,hold on
% title('Eigenvalues of R')

figure
plot(Vc(:,end-2:end)),grid on,hold on
title('Several eigenvectors estimated by eig')

%% Theoretical eigenvalues and vectors
NR = size(R,1)
d = zeros(NR,1);

SF = 1/(sqrt(NR))
V = al_gen_dftmtx(NR);
V = V*SF;

% R*V-V*D
D = V'*R*V;

figure
imagesc(real(V))
title('real of eigenvectors matrix V')

figure
imagesc(abs(D))
title('abs of eigenvalues matrix D')
% figure
% imagesc(abs(V'*V))

d = diag(D);

figure
plot(fftshift(10*log10(abs(d))),'b- o'),grid on,hold on
title('log scaled abs of eigenvalues(PSD)')
% return

figure
plot(real(V(:,1:3)))
title('real part of several eigenvectors')

figure
plot(lambda_v,'b- o'),grid on,hold on
plot(flipud(sort(abs(d))),'r- x'),grid on,hold on
title('Abs of eigenvalues of R')
return
