% Last modified on 3.14.17 by MinJu 

function a = space_avg9_mod(b)
tic

[ly, lx, lt] =size(b);
a=zeros(lx,ly,lt);

% Corners
for x=1:lx; 
    for y=1:ly;    
        if x==1 && y==1 % top-left corner 
            a(x,y,:)=(b(x,y,:)+b(x,y+1,:)+b(x+1,y,:)+b(x+1,y+1,:))/4;
        elseif x==1 && y==ly % bottom-left corner
            a(x,y,:)=(b(x,y-1,:)+b(x,y,:)+b(x+1,y-1,:)+b(x+1,y,:))/4;
        elseif x==lx && y==1 % top-right corner
            a(x,y,:)=(b(x-1,y,:)+b(x-1,y+1,:)+b(x,y,:)+b(x,y+1,:))/4;
        elseif x==lx && y==ly % bottom-right corner 
            a(x,y,:)=(b(x,y,:)+b(x,y-1,:)+b(x-1,y,:)+b(x-1,y-1,:))/4;
        end
    end
end

% Edges - corners
for y=2:ly-1; 
    for x=1:lx; 
        if x==1 % left edge 
            a(x,y,:)=(b(x,y-1,:)+b(x,y,:)+b(x,y+1,:)+b(x+1,y-1,:)+b(x+1,y,:)+b(x+1,y+1,:))/6;
        elseif x==lx % right edge 
            a(x,y,:)=(b(x-1,y-1,:)+b(x-1,y,:)+b(x-1,y+1,:)+b(x,y-1,:)+b(x,y,:)+b(x,y+1,:))/6;
        end 
    end
end

for x=2:lx-1;
    for y=1:ly; 
        if y==1 % top edge 
            a(x,y,:)=(b(x-1,y,:)+b(x-1,y+1,:)+b(x,y,:)+b(x,y+1,:)+b(x+1,y,:)+b(x+1,y+1,:))/6;
        elseif y==ly % bottom edge 
            a(x,y,:)=(b(x-1,y-1,:)+b(x-1,y,:)+b(x,y-1,:)+b(x,y,:)+b(x+1,y-1,:)+b(x+1,y,:))/6;
        end
    end
end

% Centers
for y=2:ly-1;
    for x=2:lx-1;
        a(x,y,:)=(b(x-1,y-1,:)+b(x-1,y,:)+b(x-1,y+1,:)+b(x,y-1,:)+b(x,y,:)+b(x,y+1,:)+b(x+1,y-1,:)+b(x+1,y,:)+b(x+1,y+1,:))/9;
    end 
end
toc