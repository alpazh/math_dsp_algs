close all
clear
clc

n = 60;

% Test SSVD on 1d image restoration problem
% [A,b,x] = shaw(n);
% Test SSVD on deriv2 problem
example = 3;
[A,b,x] = deriv2(n,example);

figure
plot(b),grid on
title('b')

figure
plot(x),grid on
title('x')

figure
imagesc(A),grid on
title('A')

figure
meshc(A),grid on
title('A')

% Check condition number of A,
% A is extremely ill-conditioned.
cond_A = cond(A);

% compute thin svd
[U,s,V] = csvd(A);
% define vector k values for which to compute truncated solutions
k_max = n;
k_v = (1:k_max);

%Add noise
b_exact = b;
nu = 5e-7;
b = b_exact + randn(size(b))*norm(b_exact)*nu;

% compute truncated svd
[X,rho,eta_tsvd] = tsvd(U,s,V,b,k_v);
% compute selective svd
% tau = 1e-4;  % for shaw
tau = 1e-7;% for deriv2
[x_tau,eta_ssvd,beta] = al_ssvd(U,s,V,b,tau);
err_tau = x_tau - x;

% plot absolute values of betas to compare with threshold tau
beta_abs = abs(beta);
figure
plot(log10(beta_abs),'b- o'),grid on,hold on
title('log10(abs(u''*b))')

% plot truncated solutions, compute error norms
x_norm_v = zeros(size(k_v));
err_norm_v = zeros(size(k_v));
figure(20)
plot(x,'g-'),grid on,hold on
plot(x_tau,'k-s'),grid on,hold on
figure(21)
plot(x,'g-'),grid on,hold on
plot(x_tau,'k-s'),grid on,hold on
for k = k_v%(1:4)
    xk = X(:,k);
    err_k = xk - x;
    norm_err_k = norm(err_k);
    err_norm_v(k) = norm_err_k;
    x_norm_v(k) = norm(xk);
    disp([k norm_err_k])
    figure(21)
    plot(xk,'-o'),grid on,hold on
    pause(1/24)
end

figure
plot(k_v,err_norm_v,'b- o'),grid on,hold on
plot(k_v,x_norm_v,'r- o'),grid on,hold on
legend('TSVD error norm','x{_k} norm')
xlabel('k')
ylabel('error norm and x norm for k-truncated solution')

figure
plot(k_v,err_norm_v,'b- o'),grid on,hold on
plot(k_v,ones(size(err_norm_v))*err_tau,'b- s'),grid on,hold on
legend('TSVD error norm','SSVD error norm')
xlabel('k')
ylabel('TSVD and SSVD error norm')

% figure
% loglog(x_norm_v,err_norm_v,'b- o'),grid on,hold on
% xlabel('x norm for k-truncated solution')
% ylabel('error norm for k-truncated solution')
