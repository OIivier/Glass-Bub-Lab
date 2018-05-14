function [a,b]=isoIBI_plot(N,thresh,thresh2,start,fin)
 %this function plots IBIs of all (decent) pixels over a period of time
 %(preferably 2 cycles) defined by start and fin. 
 %thresh is the threshold for activation detection, while thresh2 is the
 %threshold for gluing abberrantly short IBIs together (a.k.a IBIdoctoring)
 
sr=40;
m=size(N);
m=m(1);
a=zeros(m,m,3);
%prealloacte 3 entries for k because we only need to catch 2 successive
%IBIs and the third is in case the IBIdoctor needs to be involved.
b=zeros(m,m);
k=ones(m,m);


for i=1:m
    for j=1:m
        for t=start:fin  
           if (N(i,j,t)<thresh && N(i,j,t+1)>thresh)
                a(i,j,k(i,j))=(t/sr)+(thresh-N(i,j,t))/(sr*(N(i,j,t+1)-N(i,j,t)));
                k(i,j)=k(i,j)+1;
                if k(i,j)==3
                    b(i,j)=a(i,j,2)-a(i,j,1);
                    if b(i,j)<thresh2
                        k(i,j)=k(i,j)-1;
                        b(i,j)=0;
                    end
                end           
           end
        end
    end
end

for j=1:m
    for i=1:m
        if(b(i,j)==0)
            b(i,j)=NaN;
        end
    end
end



starttime=start/sr;
endtime=fin/sr;

 pcolor(b)
 shading flat
 h=colormap(jet);
 rgb2hsv(h);
 axis ij
 bmin=min(nonzeros(b));
 bmax=max(max(b));
 %caxis([bmin bmax])
 caxis([1 1.7])
 title(['IBIs in space from t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
 colorbar
 
 
 
 