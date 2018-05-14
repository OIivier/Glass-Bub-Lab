function [a]=isoactivation_plot(N,thresh,start,fin)
 
sr=40;
m=size(N);
m=m(1);
a=zeros(m,m);

for t=start:fin
    for j=1:m
        for i=1:m  
           if (N(i,j,t)<thresh && N(i,j,t+1)>thresh)
                a(i,j)=(t/sr)+(thresh-N(i,j,t))/(sr*(N(i,j,t+1)-N(i,j,t)));
           end
        end
    end
end

for j=1:m
    for i=1:m
        if(a(i,j)==0)
            a(i,j)=NaN;
        end
    end
end

 pcolor(a)
 shading flat
 h=colormap(jet);
 rgb2hsv(h);
 axis ij
 amin=min(nonzeros(a));
 amax=max(max(a));
 caxis([amin amax])
 title('Threshold crossing times in space')
 colorbar
 
 
 
 