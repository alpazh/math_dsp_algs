function [s_symb_hat,x_mean] = al_det_ml_pam_odd(x,OSF,M)
% function [s_symb_hat,x_mean] = al_det_ml_pam_odd(x,OSF,M)
% PAM-M (for odd value of M) ML Detector
% This is Min prob of error (MPOE) detector
% which is reduced to MAP detector
% which is reduced to ML Detector
%
% Input paramaters:
% x - input signal x,
% OSF - oversampling factor OSF (samples per symbol ratio),
% M - M different levels of PAM to detect.
%
% Output paramaters:
% s_symb_hat - detected symbols,
% x_mean - estimated mean values (per symbol).
% Test example of usage is in d\detection\test_det_ml_pam_odd.m
%
% Reference: Kay, Detection Theory, example 3.6

Nb = length(x)/OSF;
x_symb = reshape(x,OSF,Nb);
x_mean = mean(x_symb);

% calc and plot to double check slicer rule correctness
% x_mean_v = (-(M+1)/2:0.01:(M+1)/2)';
% slicer_rule = round(x_mean_v + (M+1)/2);
% slicer_rule(slicer_rule<1) = 1;
% slicer_rule(slicer_rule>M) = M;
% figure
% plot(x_mean_v,slicer_rule,'m-'),grid on
% title('Slicer Rule')

s_symb_hat = round(x_mean.' + (M+1)/2);
s_symb_hat(s_symb_hat<1) = 1;
s_symb_hat(s_symb_hat>M) = M;

return
