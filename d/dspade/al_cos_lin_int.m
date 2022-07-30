function [cos_x0_lint,x_1_i] = al_cos_lin_int(x0,lut)

% Interpolation
% x0 = 1*pi/3
x = lut(:,1);
x_1 = x(x < x0);
x_1_i = length(x_1);
% [x(x_1_i) x0 x(x_1_i+1)]
dy = lut(x_1_i+1,2)-lut(x_1_i,2);
dx = (x0-x(x_1_i))/(x(x_1_i+1)-x(x_1_i));
cos_x0_lint = lut(x_1_i,2)+dx*dy;
% [lut(x_1_i,2) cos_x0_lint lut(x_1_i+1,2)]
