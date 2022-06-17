function [x_tau,eta,beta] = al_ssvd(U,s,V,b,tau)
%AL_SSVD Selective SVD regularization.
%
% [x_tau,eta,beta] = asvd(U,s,V,b,k)
%
% Computes the selective SVD solution
%    x_tau = V(:,k_tau)*inv(diag(s(k_tau)))*U(:,k_tau)'*b .
% where k_tau is set of k-th components for which:
% norm(U(:,k)'*b) < tau

% The solution norm is returned in eta.


% Ref: Inverse Problems.
% Chapter 4. Computational Aspects: Regularization Methods, page 80

% Initialization.
% [n,p] = size(V);
[n,~] = size(V);
x_tau = zeros(n,1);
% rho = zeros(1,1);
beta = U'*b;
xi = beta./s; %xi means greek letter "ksai" not i-th x.
beta_abs = abs(beta);

% Compute
% Treat each k separately.
for k=1:n
    if (beta_abs(k)>tau)
        x_tau = x_tau + V(:,k)*xi(k);
%     else
%         rho = rho + norm(beta(k));
    end
end

eta = norm(x_tau);

% if (nargout > 1 & size(U,1) > p)
%     rho = sqrt(rho.^2 + norm(b - U(:,1:p)*beta)^2);
% end

% TODO:
% add residual norm rho = norm(A*x_tau-b) computation
return