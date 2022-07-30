close all
clear
clc

dx = 1e-1;
x_max =  1.2;
x_min = -x_max;
x1 = (x_min:dx:x_max)';
x2 = (x_min:dx:x_max)';

[X1,X2] = meshgrid(x1,x2);
F = X1.^2 + X2.^2

a1 = 1
a2 = 9
b = 1
c_x2 = (b - a1*x1)/a2;
c_z = ones(size(c_x2))*b;

figure
mesh(X1,X2,F), grid on, hold on
% figure
plot3(x1,c_x2,c_z)

x1_opt = (a1*b)/(a1.^2+a2.^2)
x2_opt = (a2*b)/(a1.^2+a2.^2)
z_opt = x1_opt.^2+x2_opt.^2
z_opt = (b*b)/(a1.^2+a2.^2)
plot3(x1_opt,x2_opt,z_opt,'r-o')

figure
plot(x1,c_x2),grid on, hold on
plot(x1_opt,x2_opt,'r *'),grid on, hold on
r = sqrt(x1_opt^2+x2_opt^2)
phi = (0:pi/36:2*pi)';
plot(r*cos(phi),r*sin(phi),'g-'),grid on, hold on


return
