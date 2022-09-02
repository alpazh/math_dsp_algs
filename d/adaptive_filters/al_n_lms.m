function [d_hat,e,w,dc] = al_n_lms(s,d,alpha,w,h_prim_path)

M = length(w);% number of taps
N = length(d);% length of the signal
epsi = 1e-6;% epsilon

u = zeros(M,1);% update vector (filter delay line)
e = zeros(N,1);% error signal
d_hat = zeros(N,1);% LMS filter output signal
dc = zeros(N,1);% mean-square deviation curve

% dc(1) = norm(w-h_prim_path).^2;

for k=1:N
    u(2:M,1) = u(1:M-1,1);
    u(1,1) = s(k);
    
    d_hat(k) = u'*w;
    
    e(k) = d(k) - d_hat(k);
    
    f = epsi + norm(u).^2;% factor
    w = w + (alpha/f)*u*e(k);
    dc(k) = norm(w-h_prim_path).^2;
end

end