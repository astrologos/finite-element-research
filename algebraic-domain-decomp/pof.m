clear;
close all;
clc;




%% Test 2 (Harwell-Boeing)
% load west0479.mat
% A = west0479;
% S = A * A' + speye(size(A));
% pct = 100 / numel(A);
% 
% figure()
% B = S
% spy(B)
% title('SSPSD Matrix Before CM Reordering')
% Bcm=(S(symrcm(S),symrcm(S)));
% figure()
% spy(Bcm)
% title('SSPSD Matrix After CM Reordering')

% %% Test 2 (airfoil)
% load airfoil
% % Scaling x and y
% x = pow2(x,-32); 
% y = pow2(y,-32);
% 
% % Forming the sparse adjacency matrix and making it positive definite
% n = max(max(i),max(j));
% A = sparse(i,j,-1,n,n);
% A = A + A';
% d = abs(sum(A)) + 1;
% A = A + diag(sparse(d));
% spy(A)
% title('SSPSD Matrix Before CM Reordering')
% B=A;
% Bcm=(A(symrcm(A),symrcm(A)));
% figure()
% spy(Bcm)
% title('SSPSD Matrix After CM Reordering')

%% Divide Into Symmetric Diagonal Blocks

n = length(Bcm)
nblocks = 10
nvert = nblocks+1;
vert = zeros([1,nvert]);
vert(1) = 1;
vert(nvert)=n;
for k=1:nvert-2
    vert(k+1) = round(((k)*n/nblocks));
end
vert;
A = Bcm;
% A1
bands = find(A(vert(2),:));
rbandwidth = bands(length(bands)) - vert(2);
A(1:vert(2)+rbandwidth,1:vert(2)+rbandwidth)=ones(vert(2)+rbandwidth);
figure(3)
spy(A)

% Ak
for k=2:nblocks-1
    bands = find(A(vert(k+1),:));
    rbandwidth = bands(length(bands)) - vert(k+1);
    A(vert(k):vert(k+1)+rbandwidth,vert(k):vert(k+1)+rbandwidth)=ones(vert(k+1)+rbandwidth-vert(k)+1);
    spy(A)
end

% AN
tempvertex=vert(length(vert)-1);
A(tempvertex:n,tempvertex:n)=ones(n-tempvertex+1);
spy(A)
title('SSPSD Matrix Divided into Diagonal Sub-Blocks after CM Reordering')