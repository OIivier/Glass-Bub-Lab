function c=movie_maker2(g)

m=length(g(1,1,:));
b=uint8(zeros(144,435,m));
c=struct('cdata', b,'colormap',[]);

for i=1:m
    pcolor(g(:,:,i));
    shading flat
    colormap(gray)
    c(i)=getframe;
end