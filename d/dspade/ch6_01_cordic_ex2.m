close all
clear
clc

phi = (0:pi/100:pi/4)';
tan_phi = tan(phi);
figure
plot(phi,tan_phi,'b- o'),grid on,hold on

y = 2.^(-7:1:0)
phi = atan(y)
phi_grad = phi*180/pi
return