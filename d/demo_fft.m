script
close all
clear all
clc
m1 = [1 2 1 2 2 2 1 1]
m2 = [1 3 1 3 3 3 1 1]
mm = [m1 m2]
emm = repmat(mm',1,1);
emm = emm(:);
emm = emm - mean(emm);
MM = fft(emm);
AMM = abs(MM);
PMM = angle(MM);

figure
plot(fftshift(AMM),'b- .'),grid on,hold on
figure
plot(fftshift(PMM),'m- .'),grid on,hold on

