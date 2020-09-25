% clc;
% clear all;
close all
addpath 'algorithms'
addpath 'graphAnalysisToolbox-1.0'
nei= 1;                % 0: 4-neighbor 1: 8-neighbor
lambda =0.1;            % smoothness parameter10
sigma_g =10;              % bw for static guidance500
sigma_u = 80;              % bw for dynamic guidance200
itr=3;%20
issparse = true;
dep = double(imread('cleandepth10.png'));
img = double(imread('cleancolor10.png'));
% I = double(imread('bigdolldepth.png'));
Iact=dep;
[m,n,d]=size(dep);
[m,n,d2]=size(img);
sigma=0.04*255;
I=dep;
I=I+sigma*randn(m,n,d);
I=I./255;
img=img./255;
[m, n,d2] = size(dep);
u0=ones(m,n,d2);
g=img;
u0=Ikmean./255;
f=u0;
[~, ~, Zg]=size(g);
[X, Y, Zu]=size(u0); N = X*Y;
[~, edges1] = lattice(X,Y,nei);
ffedges2=adjtoedges(MAP1);
edges=[ffedges2;edges1];
fVals = reshape(f,N,Zu);
if (issparse)
    A = zeros(N,1);
    A(fVals > 0) = 1;
    C=sparse(1:N,1:N,A);
    F = C*double(fVals);
else
    C = sparse(1:N,1:N,ones(N,1));
    F = double(fVals);
end
gVals = reshape(g,N,Zg);
weights_g = makeweightsL2(edges,gVals,sigma_g);%点与点的权重
gW = adjacency(edges,weights_g,N);
gW=gW-diag(diag(gW));
fprintf(1,'lambda: %d, # of steps: %d\n',lambda, itr);
WEI=WWP;
W2=WEI;
UG=edges(:,1);
VG=edges(:,2);
GD=Clusterindex(UG);
lenW=length(UG);
GW=zeros(lenW,1);
numclu=size(W2,3)
for i=1:numclu
    FF=W2(:,:,i);
    gg2=FF(VG).*(GD==i);
    GW=GW+gg2;
end
N=max(max(edges));
MAPG2=adjacency(edges,GW,N);
MAPG2=MAPG2-diag(diag(MAPG2));
N= X*Y;
for i=1:itr
    fprintf(1,'%d/%d\n',i, itr); pause(.1)
    
    %     if Zu > 1,u0 = colorspace('Lab<-', u0);end;
    uVals = reshape(u0,N,Zu);
    weights_u = makeweightsL2(edges,uVals,sigma_u);
    uW = adjacency(edges,weights_u,N);
    uW=uW-diag(diag(uW));
    %     uW=uW.*MAPG2;
    %      W=gW.*uW;
    W3=gW.*uW;
    D  = sparse(1:N,1:N,sum(W3));
    L = D-W3;
    
    R = (C+lambda*L);
    U = (R\F);
    
    u = reshape(U,X,Y,Zu);
    u0 = u;
end
figure,imshow(u0,[])
figure,imshow(Iact,[])
figure,imshow(Iact-u0.*255,[])
figure,imshow(dep,[])
Ikmean3=u0.*255;
mse = (sum(error2.^2)/(d*m*n));
maxvalue=255;
peaksnrvalue=20*log10(maxvalue)-10*log10(mse)
figure;imshow(uint8(Ikmean3),jet(255));