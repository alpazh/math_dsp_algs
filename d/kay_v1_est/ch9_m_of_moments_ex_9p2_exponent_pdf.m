script
close all
clear global
clc
% rng(123)
N = 1e4;
A = 1;
lambda = 2.7;
mu = 1/lambda
x = random('Exponential',mu,N,1);

figure
plot(x,'b- .'),grid on,hold on

figure
histogram(x)

lambda_hat = 1/mean(x)
lambda
err = lambda_hat - lambda
return
