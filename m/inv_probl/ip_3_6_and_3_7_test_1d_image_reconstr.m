close all
clear
clc

n = 24
[A,b,x] = shaw(n);

figure
plot(b),grid on
title('b')

figure
plot(x),grid on
title('x')

figure
imagesc(A),grid on
title('A')

figure
meshc(A),grid on
title('A')

% Check condition number of A,
% A is extremely ill-conditioned.
cond_A = cond(A)

[U,S,V] = svd(A);
for k = 4
    figure(1)
    plot(U(:,k),'b- o'),grid on,hold on
    % title('U vectors')
    % figure
    plot(V(:,k),'r- x'),grid on,hold off
    title('U and V vectors')
    pause(1)
end

[U,s,V] = csvd(A);
figure
eta = picard(U,s,b)
b_exact = b
b = b_exact + randn(size(b_exact))*std(b_exact)*(1e-10);
figure
eta2 = picard(U,s,b)

% Compute the partial sums of the naive SVD-based solution
b = b_exact;
[U,S,V] = svd(A);
for k = 1:n
    xk = zeros(size(x));
    for ii = 1:k
        ui = U(:,ii);
        si = S(ii,ii);
        vi = V(:,ii);
        xk = xk + ((ui'*b)/si)*vi;
    end
    err_k = xk - x;
    norm_err_k = norm(err_k);
    disp([k norm_err_k])
%     figure(100+k)
    figure(100)
    plot(x,'b- o'),grid on,hold on
    plot(xk,'r- x'),grid on,hold off
    pause(1)
    
end
figure
plot(log10(diag(S)),'b- s'),grid on
title('log10(s)')

%TODO: add 3_7 solution (noise experiment)
return
