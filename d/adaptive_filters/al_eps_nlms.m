function [y,e,w] = al_eps_nlms(s,d,mu,w)

M = length(w);% number of taps
N = length(s);% length of the signal
epsi = 1e-6;% epsilon

u = zeros(M,1);% update vector (filter delay line)
e = zeros(N,1);% error signal
y = zeros(N,1);% LMS filter output signal


for k=1:N
    u(2:M,1) = u(1:M-1,1);
    u(1,1) = s(k);
    y(k) = u.'*w;
    e(k) = d(k) - y(k);
    f = epsi + norm(u)^2;%factor
    w = w + (mu/f)*u*e(k);
end

end
