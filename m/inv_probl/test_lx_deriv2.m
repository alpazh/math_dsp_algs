close all
clear
clc
% This is the demo how noise in right hand side of Ax=b 
% affectes the decision x

% Discretize inverse problem data using regularization toolbox.
n = 3e2
example = 1%[1:3]
[A,b,x] = deriv2(n,example);
% Check condition number of A, 
% A is mildly bad conditioned.
cond_A = cond(A)

% figure
% plot(A),grid on
% title('A')

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
meshc(-A),grid on
title('A')

% Find straight forward solution x_hat and evaluate error norm.
x_hat = inv(A)*b;
err = (x_hat - x);
err_norm = norm(err)

% Add small noise in b vector and again
% find straight forward solution x_hat and evaluate error norm.
s_dist = 1/(2^14)
b_dist = b + randn(size(b))*s_dist;
x_hat_dist = inv(A)*b_dist;
err_dist = (x_hat_dist - x);
err_dist_norm = norm(err_dist)
% find residual norm
% Small residual norm does not guarantee good soltion
res_dist_norm = norm(A*x_hat_dist-b_dist)

% For the case with the noise the error is so huge that 
% the estimated x is totally useless
figure
plot(b,'b- o'),grid on,hold on
plot(b_dist,'r- x'),grid on,hold off
legend({'b','b{_dist}'})
title('True b and noised b{_dist}')

figure
plot(x_hat,'b- o'),grid on,hold on
title('Noiseless x{_hat}')
figure
plot(x_hat,'b- o'),grid on,hold on
plot(x_hat_dist,'r- x'),grid on,hold off
legend({'x_hat','x{_hat_dist}'})
title('Noiseless x{_hat} and x{_hat_dist} for the case with noise')

