function [P] = al_q_chipr2_func(nu,lambda,x,epsilon)
% function [P] = al_q_chipr2_func(nu,lambda,x,epsilon)
% computes the right-tail probability of a central or noncentral chi-squared PDF.
%
% Input paramaters:
% nu - Degree of freedom,
% lambda - Noncentrality parameter (positive), equals to 0 for central chi-squared PDF,
% x - Real scalar value of random variable,
% epsilon - max allowed error.
%
% Output paramaters:
% P - right-tail probability.
%
% Reference:
% Kay, Fundamentals of Statistical Signal Processing,
% Volume III Practical Algorithm Development

t = exp(lambda/2)*(1-epsilon);
sum_val = 1;
M = 0;
while sum_val < t
    M = M + 1;
    sum_val = sum_val + ((lambda/2)^M)/prod(1:M);
end

if rem(nu,2) == 0
    % if nu is even:
    % calculate 0-th member of the sum
    Q2 = exp(-x/2);
    g = Q2;
    for m = 4:2:nu
        g = g*x/(m-2);
        Q2 = Q2 + g;
    end
    P = exp(-lambda/2)*Q2;
    % calculate the rest of the sum
    for k = 1:M
        m = nu + 2*k;
        g = g*x/(m-2);
        Q2 = Q2 + g;
        arg_val = (exp(-lambda/2)*(lambda/2)^k)/prod(1:k);
        P = P + arg_val*Q2;
    end
else
    % if nu is odd:
    % calculate 0-th member of the sum
    P = 2*al_q_func(sqrt(x));
    % start recursion with Qchi2p_3(x)
    Q2p = sqrt(2*x/pi)*exp(-x/2);
    g = Q2p;
    if nu > 1
        for m = 5:2:nu
            g = g*x/(m-2);
            Q2p = Q2p + g;
        end
        P = P + exp(-lambda/2)*Q2p;
        % calculate the rest of the sum
        for k = 1:M
            m = nu + 2*k;
            g = g*x/(m-2);
            Q2p = Q2p + g;
            arg_val = (exp(-lambda/2)*(lambda/2)^k)/prod(1:k);
            P = P + arg_val*Q2p;
        end
    else
        % for nu == 1 0th term equals to  Qchi2_1(x) = 2*Q(sqrt(x)).
        % Add 0th and 1st terms.
        P = P + exp(-lambda/2)*(lambda/2)*Q2p;
        % calculate the rest of the sum
        for k = 2:M
            m = nu + 2*k;
            g = g*x/(m-2);
            Q2p = Q2p + g;
            arg_val = (exp(-lambda/2)*(lambda/2)^k)/prod(1:k);
            P = P + arg_val*Q2p;
        end
    end
end

return
