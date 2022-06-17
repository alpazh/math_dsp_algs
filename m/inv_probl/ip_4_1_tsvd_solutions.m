close all
clear
clc

n = 60
[A,b,x] = shaw(n);

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
cond_A = cond(A)

% compute thin svd
[U,s,V] = csvd(A);
% define vector k values for which to compute truncated solutions
k_max = n;
k_v = (1:k_max);
% compute truncated svd
[X,rho,eta] = tsvd(U,s,V,b,k_v);

% plot truncated solutions, compute error norms
x_norm_v = zeros(size(k_v));
err_norm_v = zeros(size(k_v));
figure(21)
plot(x,'g-'),grid on,hold on
for k = k_v
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
legend('error norm','x{_k} norm')
xlabel('k')
ylabel('error norm and x norm for k-truncated solution')

figure
loglog(x_norm_v,err_norm_v,'b- o'),grid on,hold on
xlabel('x norm for k-truncated solution')
ylabel('error norm for k-truncated solution')
