%script
close all
clear 
clc
% Reference:
% Detection Theory, Kay
% Chapter 2, p.48

% set x range
x = (-5:0.01:5)';
% get Q values for given x and some non Gaussian PDF
Qx = 0.5*al_q_func(x) + 0.5*al_q_func(x./sqrt(2));
% define y providing straight line for Gaussian PDF
y = [1e-4 1e-3 1e-2 1e-1 0.3 0.5 0.7 0.9 0.99 0.999 0.9999]';
% find max value for x
xmax = al_q_inv_func(min(y));
% find slope m
m = (2*min(y)-1)/(2*xmax);
% find Y values for provided slope and Q values
Y = m*al_q_inv_func(Qx)+0.5;
% define x and y for normal PDF
x_normal = (-xmax:0.01:xmax)';
y_normal = m*x_normal+0.5;

% plot and compare
figure
plot(x_normal,y_normal,'--'),grid on,hold on
plot(x,Y,'-'),grid on,hold on
xlabel('x')
ylabel('Right Tail Probability')
legend({'Q(x) (Gaussian)','non Gaussian'})
return
