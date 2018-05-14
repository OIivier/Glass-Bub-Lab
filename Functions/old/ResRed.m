function [M, Lx, Ly] = ResRed(m,dx)
% M = ResRed(m,dx), Resolution Reduction of m by taking out dx rows and columns
% from the center element. L is a cell array with the coordinates of
% elements that are being used.
[nx,ny,~] = size(m);

cx = round(nx/2);
cy = round(ny/2);

vx1 = cx : -(dx+1) : 1;
vx2 = cx : dx+1 : nx;
vx = [vx1 vx2];
vx = sort(unique(vx));

vy1 = cy : -(dx+1) : 1;
vy2 = cy : dx+1 : ny;
vy = [vy1 vy2];
vy = sort(unique(vy));

% reduced matrix
M = m(vy,vx,:);

% label matrix
% l = zeros(nx,ny);
% l(vy,vx) = 1;

% % x and y coordinates
% [row,col] = find(L==1);
% x = col;
% y = row;
% xy = [x ; y];

% coordinate matrices
[X,Y] = meshgrid(1:nx,1:ny);
Lx = X(vy,vx);
Ly = Y(vy,vx);