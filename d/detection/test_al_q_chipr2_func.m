close all
clear
clc

nu=1;lambda=2;x=0.5;epsilon=0.0001;
P = al_q_chipr2_func(nu,lambda,x,epsilon)
Pk = Qchipr2(nu,lambda,x,epsilon)
P_expected = 0.7772;
err1 = P_expected - P

nu=5;lambda=6;x=10;epsilon=0.0001;
P = al_q_chipr2_func(nu,lambda,x,epsilon);
P_expected = 0.5063;
err2 = P_expected - P

nu=8;lambda=10;x=15;epsilon=0.0001;
P = al_q_chipr2_func(nu,lambda,x,epsilon);
P_expected = 0.6161;
err3 = P_expected - P
