script
close all
clear global
clc
% ConstantScalars and Matrices
% vectors and scalar indices


%% Set parameters
% Noise deviation
Sigma = 0.05;
% Problem dimensions
M = 211;%257;
N = 211;%257;
% RNG state
rng('default')
rng(20062022)

% Time vector
t = linspace(-5,100,N);

%% Generate instrument impulse response (IIR)
g = zeros(size(1));
T0 = 10;
g(t>0) = t(t>0).*exp(-t(t>0)/T0);

Max_g = max(g);
g = g/Max_g;


%% Form matrix G
row_G_1 = zeros(1,N-1);
col_G_1 = zeros(N-1,1);
for kr=2:M
    tp=t(1)-t(kr);
    if (tp > 0)
        col_G_1(kr-1,1)=0;
    else
        col_G_1(kr-1,1)=-tp*exp(tp/T0);
    end
end
row_G_1(1,1) = col_G_1(1,1);
G = toeplitz(col_G_1,row_G_1);
% size(G)
% figure
% plot(col_G_1,'g- x'),grid on,hold off

%normalize matrix G
Delta_t = t(2)-t(1);
% G = G./(Max_g*Delta_t);
G = G./Max_g*Delta_t;

figure
colormap(jet)
imagesc(G)
H=colorbar;
set(H,'FontSize',18);
xlabel('j')
ylabel('i')
title('Matrix G')

figure
surfc(G)

%% Compute SVD of G
[U,S,V] = svd(G);

%% Generate data
% Signal magnitude
SigMagn = 2;
ScaleFactor = 1/(SigMagn^2*2);
t1 = t(1,1:end-1);
% size(t)
% size(t1)
% return
m_exact = exp(-(t1-8).^2*ScaleFactor) + 0.5*exp(-(t1-25).^2*ScaleFactor);
m_exact = m_exact.'/max(m_exact);

% Calculate ideal observed data without noise
% [size(G) size(m_exact)]
d = G*m_exact;

% Generate noised data
rng(20062022);
dn = G*m_exact+Sigma*randn(M-1,1);

figure
plot(d,'b- s'),grid on,hold on
plot(dn,'k- s'),grid on,hold on
title('Noise free data d and noised data dn')
legend({'d','dn'})

% Use whole SVD
k_trunc = N-1;
% Find Up, Vp, Sp
Up = U(:,1:k_trunc);
Vp = V(:,1:k_trunc);
Sp = S(1:k_trunc,1:k_trunc);

% Inverse solutions for noise free data d and noised data dn
m_hat_n = Vp*inv(Sp)*Up'*dn;
m_hat = Vp*inv(Sp)*Up'*d;

% size(t)
% size(g)
figure
plot(t,g,'k-s'),grid on
xlabel('Time (s)')
ylabel('IIR g')

% Singular values
figure
semilogy(diag(S),'k-o'),grid on
axis tight
xlabel('i')
ylabel('s_i')
title('Singular values')

% Solutions and true signal
figure
plot(m_exact,'b- s'),grid on,hold on
plot(m_hat,'k- o'),grid on,hold on
plot(m_hat_n,'r- x'),grid on,hold on
legend({'m_{exact}','m_{hat}','m_{hat n}'})

% Truncate SVD up to
k_trunc = 26;
% Find Up, Vp, Sp
Up = U(:,1:k_trunc);
Vp = V(:,1:k_trunc);
Sp = S(1:k_trunc,1:k_trunc);
% Inverse solutions for noise free data d and noised data dn
m_hat_tr_n = Vp*inv(Sp)*Up'*dn;
m_hat_tr = Vp*inv(Sp)*Up'*d;

% Solutions and true signal
figure
plot(m_exact,'b- s'),grid on,hold on
plot(m_hat_tr,'k- o'),grid on,hold on
plot(m_hat_tr_n,'r- x'),grid on,hold on
legend({'m_{exact}','m_{hat tr}','m_{hat tr n}'})

%% Plot L curve
m = zeros(N-1,1);
res_norm = zeros(N-1,1);
m_norm = zeros(N-1,1);

for k_tr = 1:N-1
    m = m+(U(:,k_tr)'*dn/ss(k_tr))*V(:,k_tr);
    res_norm(k_tr) = norm(G*m-dn);
    m_norm(k_tr) = norm(m);
    figure(21)
    plot(m_exact,'b- s'),grid on,hold on
    plot(m,'r-x'),grid on,hold off
    title(['k_{tr} = ' num2str(k_tr)])
    pause(1/8)
end

figure
plot(res_norm,m_norm,'b- .'),grid on,hold on
plot(res_norm(26),m_norm(26),'ro')
xlabel('Residual Norm ||Gm-d||_2')
ylabel('Solution Norm ||m||_2')
title('L curve')