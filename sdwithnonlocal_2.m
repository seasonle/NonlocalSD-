% input2 = padarray(W,[10 10],'symmetric');
% WEI=WG2;
close all
WEI=WG1;
W2=WEI;
% ff=W2(:,:,1);
[m,n,d2]=size(I);
t=10;
SS=20;
IDXm=0;
B=zeros(m*n*SS,5); 
tic
for i=1:m
     waitbar(i/m)
    for j=1:n

rmin = max(i-t,1);
rmax = min(i+t,m);
smin = max(j-t,1);
smax = min(j+t,n);

CC=(j-1)*m+i;
dd=Clusterindex(CC);
A=W2(rmin:1:rmax,smin:1:smax,dd);
A2=A'; 
[x2,x1]=meshgrid(rmin:1:rmax,smin:1:smax);
D=[x2(:),x1(:),A2(:)];
D=sortrows(D,3,'descend');
SSD=min(size(D,1),SS);
D=D(1:SSD,:);
B(IDXm+1:IDXm+SSD,:)=[i*ones(SSD,1),j*ones(SSD,1), D];  
IDXm=IDXm+SSD;
 
        
    end    
end


toc
B=B([1:IDXm],:);
U=(B(:,2)-1)*m+B(:,1);
V=(B(:,4)-1)*m+B(:,3);
nlmedge=[U,V];

SW=B(:,5);
edgesP=[U,V];
N=m*n;
% N=max(max(edgesP));

% MAP0=adjacency(edgesP,SW,N);
MAP0=sparse(U,V,SW,N,N);


MAP1=sqrt(MAP0 .* MAP0');
MAP=MAP0.*(MAP1>0);
% MAP=MAP1;



I3=I2(:,:,1);
output0=MAP*I3(:);
MAP2=sum(MAP,2);
output1=output0./MAP2;
output=reshape(output1,[m,n]);
% figure,imshow(output,[])
figure,imshow(uint8(output.*255),jet(255))

MAPS=MAP1>0;
MAPS=double(MAPS);
MAPS=sum(MAPS,2);
RR=full(MAPS);
FRR=reshape(RR,[m,n]);
figure,imshow(FRR,[]);
colormap(hot)
figure,imshow(uint8(output.*255));
colormap(hot)
figure;imshow(Ikmean,[]);colormap('hot');
figure,imshow(FRR,[]);

