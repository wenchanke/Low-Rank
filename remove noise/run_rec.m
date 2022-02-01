
clear;clc;

Image = imread('test.png') ;
Y = double(Image) ;
[M,N] = size(Y);
lambda =15; 
mu=1;
    
[X,S] = bsca_rec(Y,lambda,mu) ;

figure(1)
set(gca,'FontSize',12);
subplot(1,2,1)
imshow(Image)
xlabel(' (a) original picture')
subplot(1,2,2)
imshow(uint8(X))
xlabel('(b) noise removed')
