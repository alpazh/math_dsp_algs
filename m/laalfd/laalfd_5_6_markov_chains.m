close all
clear
clc

% Markov matrix P describes probs of transitions b/w states
P = [0.8 0.3;
     0.2 0.7]
 
% Initial state vector y0
y0 = [100 0]'

N = 30;

y = zeros(2,N)
y(:,1) = y0;
Pk = eye(2);
for k=2:N
    y(:,k) = P*y(:,k-1);
    Pk = P*Pk;
end
figure
plot(y.')
Pk
y_inf = y(:,N)

[V,D] = eig(P)
v1 = V(:,1)
v1 = v1/sum(v1)


return
