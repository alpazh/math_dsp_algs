close all
clear
clc

%% read image
filename = './images/tiger_face.jpeg';
IM = imread(filename);
IM = im2gray(IM);
figure
image(IM),colormap('gray')

%% noise
d = 0.1;
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

%% Test linearity
N = 256;
A = randn(N,N);
B = randn(N,N);
C = A+B;
AG = imguidedfilter(A);
BG = imguidedfilter(B);
CG = imguidedfilter(C);
D = C-CG;
err_norm = norm(D);
disp(['err_norm = ' num2str(err_norm)])
% figure
% imagesc(D)
return
