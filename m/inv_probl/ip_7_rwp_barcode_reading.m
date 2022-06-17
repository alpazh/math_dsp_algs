close all
clear
clc

rng(41997);
n = 256;
ksai = 1e-2
[A,b,x] = al_gen_bar_code_data(n,ksai);

figure
plot(b),grid on
title('b')

figure
plot(x,'b- o'),grid on,hold on
plot(b/max(b),'r- s'),grid on
title('x and b/max(b)')

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

figure
plot(s,'b- s'),grid on
title('Singular values')
figure
plot(log10(s),'b- s'),grid on
title('Log10 of singular values')

% %review singular vectors
% for k=1:n
%     figure(100)
%     plot(V(:,k)),grid on%plot singular vectors
%     % plot((abs(fft(V(:,k))))),grid on%plot spectrum of singular vectors
%     title(num2str(k))
%     pause(1/8)
% end

%% solve problem using regu tools

echo on

figure

% Part 1.  The discrete Picard condition
% --------------------------------------
%
% First generate a "pure" test problem where only rounding
% errors are present.  Then generate another "noisy" test
% problem by adding white noise to the right-hand side.
%
% Next compute the SVD of the coefficient matrix A.
%
% Finally, check the Picard condition for both test problems
% graphically.  Notice that for both problems the condition is
% indeed satisfied for the coefficients corresponding to the
% larger singular values, while the noise eventually starts to
% dominate.

% [A,b_bar,x] = shaw(32);
b_bar = b;
rng(41997);
e = 1e-4*randn(size(b_bar)); b = b_bar + e;
[U,s,V] = csvd(A);
subplot(2,1,1); picard(U,s,b_bar);
subplot(2,1,2); picard(U,s,b);
pause, figure



% Part 2.  Filter factors
% -----------------------
%
% Compute regularized solutions to the "noisy" problem from Part 1
% by means of Tikhonov's method and LSQR without reorthogonalization.
% Also, compute the corresponding filter factors.
%
% A surface (or mesh) plot of the solutions clearly shows their dependence
% on the regularization parameter (lambda or the iteration number).

lambda = [1,3e-1,1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5];
X_tikh = tikhonov(U,s,V,b,lambda);
F_tikh = fil_fac(s,lambda);
iter = 30; reorth = 0;
[X_lsqr,rho,eta,F_lsqr] = lsqr_b(A,b,iter,reorth,s);
subplot(2,2,1); surf(X_tikh), title('Tikhonov solutions'), axis('ij')
subplot(2,2,2); imagesc(log10(F_tikh)), axis('ij')
title('Tikh filter factors, log scale')
subplot(2,2,3); surf(X_lsqr(:,1:17)), title('LSQR solutions'), axis('ij')
subplot(2,2,4); imagesc(log10(F_lsqr(:,1:17))), axis('ij')
title('LSQR filter factors, log scale')
pause, figure

% Part 3.  The L-curve
% --------------------
%
% Plot the L-curves for Tikhonov regularization and for
% LSQR for the "noisy" test problem from Part 1.
%
% Notice the similarity between the two L-curves and thus,
% in turn, by the two methods.

subplot(1,2,1); l_curve(U,s,b); axis([1e-3,1,1,1e3])
subplot(1,2,2); plot_lc(rho,eta,'o'); axis([1e-3,1,1,1e3])
pause, figure

% Part 4.  Regularization parameters
% ----------------------------------
%
% Use the L-curve criterion and GCV to determine the regularization
% parameters for Tikhonov regularization and truncated SVD.
%
% Then compute the relative errors for the four solutions.

lambda_l = l_curve(U,s,b),   axis([1e-3,1,1,1e3]),      pause, figure
k_l = l_curve(U,s,b,'tsvd'), axis([1e-3,1,1,1e3]),      pause, figure
lambda_gcv = gcv(U,s,b),     axis([1e-6,1,1e-9,1e-1]),  pause, figure
k_gcv = gcv(U,s,b,'tsvd'),   axis([0,20,1e-9,1e-1]),    pause, figure

x_tikh_l   = tikhonov(U,s,V,b,lambda_l);
x_tikh_gcv = tikhonov(U,s,V,b,lambda_gcv);
if isnan(k_l)
    x_tsvd_l = zeros(32,1); % Spline Toolbox not available.
else
    x_tsvd_l = tsvd(U,s,V,b,k_l);
end
x_tsvd_gcv = tsvd(U,s,V,b,k_gcv);

% figure
plot(x,'b- o'),grid on,hold on
plot(x_tsvd_l,'r- s'),grid on,hold on
plot(x_tsvd_gcv,'m- d'),grid on,hold on
title('TSVD solutions')

figure
plot(x,'b- o'),grid on,hold on
plot(x_tikh_l,'r- s'),grid on,hold on
plot(x_tikh_gcv,'m- d'),grid on,hold on
title('Tikhonov solutions for lambda_l and lambda_gcv')

norm(x-x_tikh_l)/norm(x)
norm(x-x_tikh_gcv)/norm(x)
norm(x-x_tsvd_l)/norm(x)
norm(x-x_tsvd_gcv)/norm(x)

% disp([norm(x-x_tikh_l),norm(x-x_tikh_gcv),...
%  norm(x-x_tsvd_l),norm(x-x_tsvd_gcv)]/norm(x))

lambda = 10.^(-5:0.125/2:0)%
X_tikh = tikhonov(U,s,V,b,lambda);
for k = 1:size(X_tikh,2)
    figure(201)
    plot(x,'b- o'),grid on,hold on
    plot(X_tikh(:,k)),grid on,hold off
    title(['Tikhonov solution for lambda = ' num2str(lambda(k))])
    pause(1/2)
end

return
