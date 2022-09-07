function [thresh] = al_compute_energy_detector_thresh(N,Pfa,var_wgn,x_max)
% [thresh] = al_compute_energy_detector_thresh(N,Pfa,var_wgn,x_max)
% Computes threshold gamma for Energy Detector
%
% Input paramaters:
% N - length of input signal,
% var_wgn - Noise Variance,
% Pfa - required Probability of the False Alarm,
% x_max - maximum value of the signal.
%
% Output paramaters:
% thresh - Threshold gamma for the given Pfa.
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development, equation (10.13)

x_min = 0;
err = 1000;
while abs(err) > Pfa/100
    x = (x_min + x_max)/2;
    Pfa_hat = al_q_chipr2_func(N,0,x,Pfa/100);
    err = Pfa_hat - Pfa;
    if err < 0
        x_max = x;
    else
        x_min = x;
    end
end
thresh = var_wgn*x;
return
