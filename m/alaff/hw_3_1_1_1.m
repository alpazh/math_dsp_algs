close all
clear
clc

% Create Vandermonde matrix for m points at interval [x_min x_max]
% for order p from 0 to Nmax and estimate condition number
Nmax = 4;
m = 5000;
x_min = 0;
x_max = 1;
cond_p_V = zeros(Nmax,1);
cond_p_L = zeros(Nmax,1);
for p=1:Nmax
    [V] = create_vandermonde(m,p,x_min,x_max);
    [L] = create_legendre(m,p,x_min,x_max);
    cond_p_V(p) = cond(V);
    cond_p_L(p) = cond(L);
end
% Plot condition number vs polynomial order p
figure
plot((1:Nmax)-1,log10(cond_p_V),'r- x'),grid on,hold on
plot((1:Nmax)-1,log10(cond_p_L),'b- o'),grid on
title('Condition number of Vandermonde matrix and Legendre matrix')
legend('Vandermonde','Legendre')

return
