close all
clear
clc

n = 32;
example = 1;
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

figure
plot(s,'b- s'),grid on
title('Singular values')
figure
plot(log10(s),'b- s'),grid on
title('Log10 of singular values')
% Add a small amount of noise to the right-hand side
b_exact = b;
nu_e = 1e-3;
e = nu_e*randn(size(b));
b = b_exact + e;

lambda_v = logspace(1,-5,20)
% for kl = 1:length(lambda_v)
%     lambda = lambda_v(k)
X = tikhonov(U,s,V,b,lambda_v);
% end

method = 'Tikh';
f = fil_fac(s,lambda_v,method);

figure
surfc(f),grid on
title('Filter factors')
% legend('1e-3','1e-2','1e-1','1')

figure
surfc(X),grid on
title('Tikhonov Solutions')
% legend('1e-3','1e-2','1e-1','1')

nl = size(X,2);
norm_err = zeros(nl,1);
norm_x_v = zeros(nl,1);
norm_res_v = zeros(nl,1);
for kl = 1:nl
    x_tau_k = X(:,kl);
    norm_err(kl) = norm(x_tau_k-x);
    norm_x_v(kl) = norm(x_tau_k);
    norm_res_v(kl) = norm(A*x_tau_k-b);
end
figure
plot(lambda_v,norm_err,'b- s'),grid on
title('Error norm')
% figure
% plot(lambda_v,log10(norm_x_v),'b- s'),grid on
% title('norm x')
figure
plot(lambda_v,log10(norm_res_v),'b- s'),grid on
title('norm Ax-b')

figure
loglog(norm_res_v,norm_x_v,'b- o'),grid on
ylabel('||x||')
xlabel('||Ax-b||')
title('L curve')

figure
[reg_corner,rho,eta,reg_param] = l_curve(U,s,b,method);
return
