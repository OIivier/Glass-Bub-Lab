function a=space_avg9(b)
tic
m=size(b);
ly=m(1);
lx=m(2);
lt=m(3);
a=zeros(lx-2,ly-2,lt);

for y=2:ly-1;
    for x=2:lx-1; 
        a(x,y,:)=(b(x-1,y-1,:)+b(x-1,y,:)+b(x-1,y+1,:)+b(x,y-1,:)+b(x,y,:)+b(x,y+1,:)+b(x+1,y-1,:)+b(x+1,y,:)+b(x+1,y+1,:))/9;
    end 
end
toc