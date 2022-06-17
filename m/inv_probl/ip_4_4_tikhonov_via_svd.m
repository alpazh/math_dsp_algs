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
e = 1e-3*randn(size(b));
b = b_exact + e;

lambda_v = 10.^(-3:1:0);
% for kl = 1:length(lambda_v)
%     lambda = lambda_v(k)
X = tikhonov(U,s,V,b,lambda_v);
% end

method = 'Tikh';
f = fil_fac(s,lambda_v,method);

figure
plot(f,'-v'),grid on
title('Filter factors')
legend('1e-3','1e-2','1e-1','1')

figure
plot(X,'-s'),grid on
title('Tikhonov Solutions')
legend('1e-3','1e-2','1e-1','1')
return
