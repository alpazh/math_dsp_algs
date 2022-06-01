function [X] = create_legendre(m,n,x_min,x_max)
% create_vandermonde matrix for given size m, order n on interval
% from x_min to x_max
dx = (x_max-x_min)/(m-1);
x = (x_min:dx:x_max)';
X = zeros(m,n);
for k=1:n
    switch k
        case 1
            v = ones(m,1);
        case 2
            v = 2*x+1;
        case 3
            v = 6*x.^2-6*x+1;
        case 4
            v = 20*x.^3-30*x.^2+12*x-1;
    end
    X(:,k) = v;
end
end

