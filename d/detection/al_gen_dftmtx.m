function F = al_gen_dftmtx(N)

F = zeros(N,N);
tv = (0:N-1)';
for i=0:N-1
    fi = i/N;
    F(:,i+1) = exp(-1j*2*pi*fi*tv);
end

return
