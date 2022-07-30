function [x_gmm] = al_gen_gmm(N,epsln,sigma1,sigma2)
% Generate Gausian Mixture Model data
phi1 = randn(N,1)*sigma1;
phi2 = randn(N,1)*sigma2;
selctr = (rand(N,1) > (1-epsln)) + 1;
x_gmm = phi1;
x_gmm(selctr == 2) = phi2(selctr == 2);
return