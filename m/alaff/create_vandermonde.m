function [X] = create_vandermonde(m,n,x_min,x_max)
% create_vandermonde matrix for given size m, order n on interval 
% from x_min to x_max
dx = (x_max-x_min)/(m-1);
x = (x_min:dx:x_max)';
X = x.^(0:n);
end

