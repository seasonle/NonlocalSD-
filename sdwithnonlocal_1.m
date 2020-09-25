% clc;
% clear;
close all;
%
% 
filename1='cleancolor10.png';
filename2='cleandepth10.png';
I=double(imread(filename2));
G = imread(filename1); 
G=rgb2gray(G);
G=double(G);
img = G;

Iact=I;
[m,n,d]=size(I);
[m,n,d2]=size(G);
sigma=0.01*250;
S=10;K=3;
pcadim=6;
fast_flag=1;
Iact=I;
I=I+sigma*randn(m,n,d);
G=G+0*randn(m,n,d2);
I2=I./255;
G=G./255;

error2 = reshape(Iact-I, [d*m*n,1]);
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
max(I(:));
peaksnrvalue=20*log10(maxvalue)-10*log10(mse);
fprintf('\n The Peak-SNR value for noisy image is %0.4f (db) \n', peaksnrvalue);

[h1 ,w1 ,nc1] = size(I2);
[h2 ,w2 ,nc2] = size(G);
radius=3;
Apca=compute_pca(I2, K, pcadim);
nneighbors = (2*radius + 1)^2;%7*7

H1 = zeros([h1 w1 nneighbors*nc1]);%512*512*（7*7*3）
H2 = zeros([h2 w2 nneighbors*nc2]);%512*512*（7*7*3）
n = 1;

for i = -radius:radius
    for j = -radius:radius
        dist2  = i^2 + j^2;
        weight = exp(-dist2 / 2 / (radius / 2));%高斯核加权  为什么要这样移位加权？？？我的理解是一个点非局部(这里注释错)邻域就是这些点，点与点之间做高斯加权
%         weight=1;
        C1 = circshift(I2, [i j]);%行i 列j移位
        C2 = circshift(G, [i j]);%行i 列j移位
        H1(:,:, nc1*(n-1)+1:nc1*n) = C1 * weight;
        H2(:,:, nc2*(n-1)+1:nc2*n) = C2 * weight;
        n = n + 1;
    end
end


 H1 = reshape(H1, h1*w1, nneighbors*nc1);
    H1 = bsxfun(@minus, H1, mean(H1));%样本中心化
%  H1 = H1 - repmat(mean(H1,1), h1*w1, 1);
 H2 = reshape(H2, h2*w2, nneighbors*nc2);
%  H2 = H2 - repmat(mean(H2,1), h2*w2, 1);
    H2 = bsxfun(@minus, H2, mean(H2));%样本中心化   

X=H1;
Y=H2;
PCD=6;
[XL,YL] = plsregress(Y,X,PCD);
H1C=H1*YL;
H1C=reshape(H1C,[h1 w1 PCD]);
H2C=H2*XL;
H2C=reshape(H2C,[h2 w2 PCD]);

I1C=sum(H1C,3)./PCD;
I2C=sum(H2C,3)./PCD;
figure,imshow(I1C,[])
figure,imshow(I2C,[])
Cluster=15;
pcadim=3;
pd=6;

Apca(:,:,1:3)=H1C(:,:,1:3).*1;
Apca(:,:,4:6)=H2C(:,:,1:3).*0.01;

[m,n,d]=size(I);
pcadim=size(Apca,3);
Ares=reshape(Apca,m*n,pcadim);
[Centre,mincentre]=kmeans_recursive(Ares,Cluster);
Clusterindex=reshape(mincentre,m,n);



% Apcares=imresize(Apca,[256 256]);
pd=size(Apca,3)
% Ares=reshape(Apcares,size(Apcares,1)*size(Apcares,2),pd);
% [Centre,mincentre]=kmeans_recursive(Ares,Cluster);
% ff=reshape(mincentre,256,256);
stestpca=sum(Apca,3)./pd;
figure,imshow(stestpca,[])
% Filtering
spatialtype='box';     
convmethod='O1'; % 'matlab' for matlab convolutions and 'O1' for O(1) convolutions

[Ikmean ,W3,WG1,C1]=fastKmeansfiltapproxinter4(I2,S,150*5/256,Centre,spatialtype,convmethod,fast_flag,Apca);%depth
% [Ikmean ,W3,WG1,C1]=fastKmeansfiltapproxinter4(I2,S,2000*sigma/255,Centre,spatialtype,convmethod,fast_flag,Apca);
% [Ikmean ,W,WG,C1]=fastKmeansfiltapproxinter4(I2,S,100,Centre,spatialtype,convmethod,fast_flag,Apca);%150
Ikmean=Ikmean.*255;
Tkmeans=toc;
fprintf('non local means by Kmeans complete with %d clusters \n',Cluster);
fprintf('time for non local means (ms)=%3.0f \n',Tkmeans*1000);

%% Displaying noisy and filtered image
figure;
imshow(uint8(Iact)),%title('Original image');
figure;
imshow(uint8(I)),%title('input image');
figure;
imshow(uint8(Ikmean));%title(['fast non local means with ',num2str(Cluster),' clusters']);
figure;
G=G.*255;
imshow(uint8(G)),%title('input image');
 figure;imshow(Ikmean,[]);colormap('hot');

%% PSNR calculation
[m,n,d]=size(I);
error2 = reshape(Iact-Ikmean, [d*m*n,1]);
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
peaksnrvalue=20*log10(maxvalue)-10*log10(mse);
ssimval=ssimcalculate(Iact,Ikmean);
fprintf('\n mean sq error=%f \t The Peak-SNR value for filtered image is %0.4f (db) \t SSIM value is %0.4f \n',sqrt(mse), peaksnrvalue, ssimval);

