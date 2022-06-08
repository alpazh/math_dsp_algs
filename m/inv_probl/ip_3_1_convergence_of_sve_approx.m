close all
clear
clc

powers_of_n = [2:10];
example = 1;
for p = powers_of_n
    disp('__________________')
    p
    n = 2^p
    [A,b,x] = deriv2(n,example);
%     [A,b,x] = gravity(n);
    [U,S,V] = svd(A);
    s = diag(S);
    figure(p)
    plot(log10(s),'b- o'),grid on,hold on
    
    for j = (1:4)
        figure(100*p+j)
        h = 1/n; stairs((0:n)*h,[U(:,j);U(n,j)]/sqrt(h))
    end
end

return
