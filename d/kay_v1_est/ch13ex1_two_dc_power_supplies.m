script
close all
clear all
clc
% Kay, Estimation Theory, Chapter 13, p.427, 
% Example 13.1 Two DC Power Supplies


%% 2nd order State Model:
A0 = 5;
p = 2;
r = 2;
N = 3e2;
sigma_u = 0.01;
sigma_s1 = 0.01;
sigma_s2 = 0.01;
mu_u = 0;
mu_s = ones(p,1)*A0;

a1 = 0.999;
a2 = 0.998;
A = [a1 0; 0 a2]
B = eye(2)
Cs = [sigma_s1 0; 0 sigma_s2]
UCs = chol(Cs)
LCs = UCs.'

s_v = zeros(p,N);
s_1 = LCs*randn(p,1) + mu_s;
u_v = randn(r,N)*sigma_u + mu_u;
s_v(:,1) = s_1;
for k = 2:N
    u = u_v(:,k);
    s = A*s_1 + B*u;
    s_1 = s;
    s_v(:,k) = s;
end

figure
plot(s_v'),grid on

return
