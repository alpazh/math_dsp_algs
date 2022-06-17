function [A,b,x] = al_gen_bar_code_data(n,ksai,upf)

if nargin < 2
    ksai = 0.01;
end

if nargin < 3
    upf = 8;
end


A = zeros(n,n);

ni = 1/n;
ksai2n2 = ksai^2*n^2;

for ii=1:n
    for jj=1:n
        A(ii,jj) = ni*exp(-(ii-jj)^2/ksai2n2);
    end
end

nb = ceil(n/upf);
xi = (rand(nb,1) > 0.5)*1.0;
xi = repmat(xi,1,upf);
xi = xi';
xi = xi(:);
x = xi(1:n,1);

b = A*x;
return