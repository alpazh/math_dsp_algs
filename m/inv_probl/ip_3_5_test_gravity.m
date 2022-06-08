close all
clear
clc
% This is the demo how noise in right hand side of Ax=b 
% affectes the decision x

% Discretize inverse problem data using regularization toolbox.
n = 2^6;
example = 1;%[1:3]
a_interv = 0;
b_interv = 1;
d = 0.25/2^1;
[A,b,x] = gravity(n,example,a_interv,b_interv,d);
% Check condition number of A, 
% A is bad conditioned.
cond_A = cond(A);

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
% For this problem even double-precision quantization noise is too much
% to get any meaningful solutions
x_hat = inv(A)*b;
err = (x_hat - x);
err_norm = norm(err);

% Add small noise in b vector and again
% find straight forward solution x_hat and evaluate error norm.
s_dist = 1/(2^14);
b_dist = b + randn(size(b))*s_dist;
x_hat_dist = inv(A)*b_dist;
err_dist = (x_hat_dist - x);
err_dist_norm = norm(err_dist);
% find residual norm
% Small residual norm does not guarantee good soltion
res_dist_norm = norm(A*x_hat_dist-b_dist);

% For the case with the noise the error is so huge that 
% the estimated x is totally useless
figure
plot(b,'b- o'),grid on,hold on
plot(b_dist,'r- x'),grid on,hold off
legend({'b','b{_dist}'})
title('True b and noised b{_dist}')

figure
plot(x,'b- o'),grid on,hold on
title('x')
figure
plot(x,'b- o'),grid on,hold on
plot(x_hat,'r- x'),grid on,hold on
legend({'x','x_hat'})
title('Noiseless x{_hat}')
figure
plot(x_hat,'b- o'),grid on,hold on
plot(x_hat_dist,'r- x'),grid on,hold off
legend({'x_hat','x{_hat_dist}'})
title('Noiseless x{_hat} and x{_hat_dist} for the case with noise')

d_v = (0.05:0.025:0.3)';
cond_A_v = zeros(size(d_v));
for k = 1:length(d_v)
    d = d_v(k);
    [A,b,x] = gravity(n,example,a_interv,b_interv,d);
    % Check condition number of A, 
    % A is bad conditioned.
    cond_A_v(k) = cond(A);
end

figure
plot(d_v,log10(cond_A_v),'b- s'),grid on
xlabel('d'),ylabel('log10(cond(A))')

d = 0.15;
n_v = 2.^(2:10)';
cond_A_v = zeros(size(n_v));
for k = 1:length(n_v)
    n = n_v(k);
    [A,b,x] = gravity(n,example,a_interv,b_interv,d);
    % Check condition number of A, 
    % A is bad conditioned.
    cond_A_v(k) = cond(A);
end

figure
plot(log2(n_v),log10(cond_A_v),'b- s'),grid on
xlabel('log2(n)'),ylabel('log10(cond(A))')

d_v = (0.05:0.025:0.25)';
np = (2:10)';
n_v = 2.^np;
mat_cond_A = zeros(size(d_v,1),size(n_v,1));
for kd = 1:length(d_v)
    for kn = 1:length(n_v)
        d = d_v(kd);
        n = n_v(kn);
        [A,b,x] = gravity(n,example,a_interv,b_interv,d);
        % Check condition number of A, 
        % A is bad conditioned.
        mat_cond_A(kd,kn) = cond(A);
    end
end

figure
surf(log10(mat_cond_A))
