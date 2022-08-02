function s_hat = al_mmse_chan_est(Cs,H,sigma,x)
% Kay, Estimation Theory, Chapter 11.7, p.367, eq. 11.39:
% MMSE Signal Estimation from the output of the known channel

N = length(x);
I = eye(N);
MTI = H*Cs*H'+sigma^2*I;
% cond_HCHt = cond(H*Cs*H')
% cond_MTI = cond(MTI)
s_hat = Cs*H'*inv(MTI)*x;

return
