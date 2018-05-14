close all

x=50:70;
t=1498/40:1/40:1794/40;
y=zeros(length(x),length(t));
for i=1:length(x)
    for j=1:length(t)
        y(i,j)=N1(40,i,j);
    end
end

surface(t,x,y)
shading interp