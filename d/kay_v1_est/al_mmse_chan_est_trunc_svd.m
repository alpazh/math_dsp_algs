function s_hat = al_mmse_chan_est_trunc_svd(Cs,H,x,r_trunc)
% Kay, Estimation Theory, Chapter 11.7, p.367, eq. 11.39:
% Signal Estimation from the output of the known channel with
% truncated SVD

HCHt = (H*Cs*H');
% cond_HCHt = cond(H*Cs*H')
[U,S,V] = svd(HCHt);
sv = diag(S);
svi = 1./sv;
svi_tr = svi(1:r_trunc,1);
HCHt_i_tr = V(:,1:r_trunc)*diag(svi_tr)*U(:,1:r_trunc)';
% cond_MTI = cond(MTI)
s_hat = Cs*H'*HCHt_i_tr*x;

return
