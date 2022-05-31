close all

filename = 'panda1.jpg';
M=imread(filename);
size(M)
GS = sum(M,3)/(3*256);
filename_g = [filename(1:end-4) '_g.jpg']
imwrite(GS,filename_g);

GS=imread(filename_g);
size(GS)

GSd = double(GS);
[U,S,V] = svd(GSd);


N=42
for k=1:N
    GSk = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    figure(100+k)
    imagesc(GSk),colormap('gray')
    pause
end

figure
imagesc(GS),colormap('gray')

figure
plot(20*log10(diag(S)/S(1,1))),grid on
return
