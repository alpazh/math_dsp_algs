function y = al_q_inv_func(x)

% y = al_q_func(x) returns inverse Q function, which is
% the value y for which N(0,1) > y with probability equal to x

y = sqrt(2)*erfinv(1-2*x);

return
