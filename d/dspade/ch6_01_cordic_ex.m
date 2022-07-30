close all
clear
clc
rng(20)
z = randn(1,1)+1j*randn(1,1);

x = real(z);
y = imag(z);
A = abs(z);
phi = angle(z)*180/pi;

cos_phi = cos(angle(z));
cos_phi_ratio = x/A;

figure
plot([0 real(z)],[0 imag(z)],'k- o'),grid on,hold on

% compute cos(x)
% Look Up Table
phi_step = pi/2
x = [0:phi_step:2*pi]';
lut = [x cos(x)]

% Interpolation
x0 = 3*pi/8
[cos_x0_lint,x_1_indx] = al_cos_lin_int(x0,lut)
% [x(x_1_indx) x0 x(x_1_indx+1)]
% [lut(x_1_indx,2) cos_x0_lint lut(x_1_indx+1,2)]
cos_x0 = cos(x0)

% return
% Test LUT with lin int
phi_step = pi/32
x_t = [phi_step:phi_step:2*pi-phi_step]';

for k=1:length(x_t)
    xk = x_t(k);
    [cos_xk_lint,x_1_indx] = al_cos_lin_int(xk,lut);
    cos_x_t_lint(k) = cos_xk_lint;
    cos_x_t(k) = cos(xk);
end

figure
plot(cos_x_t_lint,'r- s'),grid on,hold on
plot(cos_x_t,'b-o'),grid on,hold on





