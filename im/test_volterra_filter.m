% close all
clear
clc

%% read image
filename = './images/tiger_face.jpeg';
IM = imread(filename);
IM = im2gray(IM);
% IM = IM(1:700/2,1:700/2);
figure
image(IM),colormap('gray')

%% noise
d = 0.5;
% d = 1.0;
IMN = imnoise(IM,'salt & pepper',d);
figure
image(IMN),colormap('gray')


%% guided filter
IMMF = medfilt2(IMN,[21 21]);
IMF = imguidedfilter(IMN,IMMF);
% IMF = imguidedfilter(IMN,G);
figure
image(IMMF),colormap('gray')
figure
image(IMF),colormap('gray')

%% Identify Volterra Filter Coefficients
S = double(IM);
N = double(IMN);
Sh = size(S,1)
Sw = size(S,2)
% return
r = 2
sz = 2*r+1
fl = sz.^2

FiltOrd = 2
neq = 0
Xk_len = zeros(FiltOrd,1);
for k=1:1:FiltOrd
    %(factorial(fl+k-1)/(factorial(k)*factorial(fl-1)))
    k
    neq_k = (factorial(fl+k-1)/(factorial(k)*factorial(fl-1)))
    neq = neq + neq_k;
    Xk_len(k) = neq_k;
end
neq

EX = zeros(neq,1);

k_ord = 1;
X1_ind = zeros(Xk_len(k_ord),k_ord);
NEL=1;
% X1_ind = zeros(Xk_len)
for I1=1:fl
%     X1(NEL)=fl-I1+1;
    X1(NEL)=I1;
    X1_ind(NEL)=fl-I1+1;
    NEL=NEL+1;
end
X1_ind

k_ord = 2;
X2_ind = zeros(Xk_len(k_ord),k_ord);
NEL=1;
for I1=1:fl
    for I2=I1:fl
        % X2(NEL)=X1(I1)*X1(I2);
        X2(NEL)=10*(I1) + (I2);
        X2_ind(NEL,:) = [I1 I2];
        NEL=NEL+1;
    end
end
X2_ind
X = [X1 X2];

% k_ord = 3;
% X3_ind = zeros(Xk_len(k_ord),k_ord);
% size_X3_ind = size(X3_ind)
% NEL=1;
% for I1=1:fl
%     for I2=I1:fl
%         for I3=I2:fl
%         % X3(NEL)=X1(I1)*X1(I2)*X1(I3);
%         X3(NEL)=100*(I1) + 10*(I2) + I3;
%         X3_ind(NEL,:) = [I1 I2 I3];
%         NEL=NEL+1;
%         end
%     end
% end
% X3_ind

% X = [X1 X2 X3];

X = (1:fl)
Xp = zeros(1,neq);
size(Xp)
% Xp = [X(X1_ind) X(X2_ind(:,1)).*X(X2_ind(:,2)) X(X3_ind(:,1)).*X(X3_ind(:,2)).*X(X3_ind(:,3))]
Xp = [X(X1_ind) X(X2_ind(:,1)).*X(X2_ind(:,2))]
size(Xp)

% k_ord = 4;
% X4_ind = zeros(Xk_len(k_ord),k_ord);
% size_X4_ind = size(X4_ind)
% NEL=1;
% for I1=1:fl
%     for I2=I1:fl
%         for I3=I2:fl
%             for I4=I3:fl
%             % X4(NEL)=X1(I1)*X1(I2)*X1(I3)*X1(I4);
%             X4(NEL)=1000*(I1) + 100*(I2) + 10*I3 + I4;
%             X4_ind(NEL,:) = [I1 I2 I3 I4];
%             NEL=NEL+1;
%             end
%         end
%     end
% end
% X4_ind
% X = [X1 X2 X3 X4];
neq
XtX=zeros(neq,neq);
B=zeros(neq,1);
Xmask = zeros(fl,fl);
X = zeros(1,neq);
X2_ind1 = X2_ind(:,1);
X2_ind2 = X2_ind(:,2);
for kh = r+1:Sh-r
    for kw = r+1:Sw-r
        yn = S(kh,kw);
        Xmask = N(kh-r:kh+r,kw-r:kw+r);
%         Xmask = N(kh-r:kh+r,kw-r:kw+r).';
        % figure(101)
        % image(Xmask),colormap('gray')
        % pause(1/24)

        Xin = (Xmask(:))';
%         X = [Xin(X1_ind) Xin(X2_ind(:,1)).*Xin(X2_ind(:,2)) Xin(X3_ind(:,1)).*Xin(X3_ind(:,2)).*Xin(X3_ind(:,3))];
%         X = [Xin(X1_ind) Xin(X2_ind(:,1)).*Xin(X2_ind(:,2))];
%         X = [Xin(X1_ind) Xin(X2_ind1).*Xin(X2_ind2)];
        X(1,1:fl) = Xin(X1_ind);
        X(1,fl+1:end) = Xin(X2_ind1).*Xin(X2_ind2);
        
        Xt = X';
        B_k=Xt*yn;
        XtX_k=Xt*X;
        XtX=XtX+XtX_k;
        B=B+B_k;
    end
end
XtX = XtX+eye(neq)*1e3;
cond_XtX = cond(XtX)
H1 = inv(XtX)*B;
H2 = pinv(XtX)*B;
H3 = XtX\B;
neq
S = svd(XtX);
figure
plot(20*log10(S),'b- .'),grid on,hold on

figure
plot((H1),'r- s'),grid on,hold on
plot((H2),'g- o'),grid on,hold on
plot((H3),'b- d'),grid on,hold on

figure
plot(H1(fl+1:end),'r- s'),grid on,hold on
plot(H2(fl+1:end),'g- o'),grid on,hold on
plot(H3(fl+1:end),'b- d'),grid on,hold on

disp('VNN filter...')
filename = './images/tiger_face.jpeg';
IM = imread(filename);
IM = im2gray(IM);



F = zeros(size(N));
for kh = r+1:Sh-r
    for kw = r+1:Sw-r
        Xmask = N(kh-r:kh+r,kw-r:kw+r);
%         Xmask = N(kh-r:kh+r,kw-r:kw+r).';

        Xin = (Xmask(:))';
%         X = [Xin(X1_ind) Xin(X2_ind(:,1)).*Xin(X2_ind(:,2)) Xin(X3_ind(:,1)).*Xin(X3_ind(:,2)).*Xin(X3_ind(:,3))];
%         X = [Xin(X1_ind) Xin(X2_ind(:,1)).*Xin(X2_ind(:,2))];
%         X = [Xin(X1_ind) Xin(X2_ind1).*Xin(X2_ind2)];
        X(1,1:fl) = Xin(X1_ind);
        X(1,fl+1:end) = Xin(X2_ind1).*Xin(X2_ind2);

        F(kh,kw) = X*H2;
    end
end

figure
image(N),colormap('gray')
figure
image(F),colormap('gray')
return
