close all
clear
clc

N = 8
Pfa = 0.1
var_wgn = 1.0
x_max = 100
thresh = al_compute_energy_detector_thresh(N,Pfa,var_wgn,x_max)
thresh_expected = 13.37890625
err1 = thresh_expected - thresh
