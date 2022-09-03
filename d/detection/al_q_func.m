function y = al_q_func(x)

% y = al_q_func(x) returns probability of the N(0,1) PDF right tail 
% for given value x (i.e. probality that N(0,1)>x).

y = 0.5*erfc(x/sqrt(2));

return
