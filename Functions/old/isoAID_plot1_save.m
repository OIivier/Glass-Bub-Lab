function [a,b,c,g]=isoAID_plot1_save(N,thresh,thresh2,start,fin)

%written by Bartek Borek in May 2009

 %this function calculates and plots activation times (=a),IBIs (=b), and change in IBIs (=c) of all (decent) pixels over a period of time
 %(preferably 3 cycles) defined by start and fin. 
 %thresh is the threshold for activation detection, while thresh2 is the
 %threshold for gluing abberrantly short IBIs together (a.k.a IBIdoctoring)
 
sr=40;
starttime=start/sr;
endtime=fin/sr;
m=size(N);
m=m(1);
a=zeros(m,m,3);
%prealloacte 3 entries for k because we only need to catch 2 successive
%IBIs.
b=zeros(m,m,2);
c=zeros(m,m);

for i=1:m
    for j=1:m
        %set counters for each pixel
        k=1; 
        l=1;
        for t=start:fin  
           if (N(i,j,t)<thresh && N(i,j,t+1)>thresh)
                a(i,j,k)=(t/sr)+(thresh-N(i,j,t))/(sr*(N(i,j,t+1)-N(i,j,t)));
                
                if k==1 %this checks if your activation time is a abnormally late, and sets it to zero (undetected) if it is. 
                    if (a(i,j,k))>(start/sr)+((fin-start)/(3*sr))
                        a(i,j,k)=0;
                    end
                end
                if k==2
                    if (a(i,j,k))>(start/sr)+(2*(fin-start)/(3*sr))
                        a(i,j,k)=0;
                    end
                end
                
                k=k+1;
                
                if k>2
                    b(i,j,l)=a(i,j,l+1)-a(i,j,l);
                    if b(i,j,l)<thresh2
                        k=k-1;
                        b(i,j,l)=0; %this way if you cut after doctor, but before full IBI it let's you know.
                    else 
                        l=l+1;
                    end
                end           
           end
        end
        
        if (b(i,j,1)==0 || b(i,j,2)==0)
            c(i,j)=0;
        else
            c(i,j)=b(i,j,2)-b(i,j,1);
        end
        
    end
end



for j=1:m
    for i=1:m
        if(c(i,j)==0)
            c(i,j)=NaN;
        end
        for s=1:3
            if(b(i,j,s)==0)
                b(i,j,s)=NaN;
            end
            if(a(i,j,s)==0)
                a(i,j,s)=NaN;
            end
        end
    end
end

a1=a(:,:,1);
%figure
%subplot(3,3,1)
pcolor(a1)
shading flat
h=colormap(jet);
rgb2hsv(h);
axis ij
%a1min=min(nonzeros(a1));
%a1max=max(max(a1));
a1min=start/sr;
a1max=a1min+((fin-start)/(3*sr));
caxis([a1min a1max])
title(['CT for the first activation during t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
colorbar
 
a2=a(:,:,2);
figure
%subplot(3,3,2)
pcolor(a2)
shading flat
h=colormap(jet);
rgb2hsv(h);
axis ij
a2min=a1max;
a2max=a2min+((fin-start)/(3*sr));
caxis([a2min a2max])
title(['CT for the second activation during t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
colorbar

a3=a(:,:,3);
figure
%subplot(3,3,3)
pcolor(a3)
shading flat
h=colormap(jet);
rgb2hsv(h);
axis ij
a3min=a2max;
a3max=a3min+((fin-start)/(3*sr));
caxis([a3min a3max])
title(['CT for the third activation during t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
colorbar


b1=b(:,:,1);
figure
%subplot(3,3,4)
pcolor(b1)
shading flat
h=colormap(jet);
rgb2hsv(h);
axis ij
%bmin=min(nonzeros(b1));
%bmax=max(max(b1));
%caxis([bmin bmax])
caxis([0.4 2])
title(['IBIs for first two activations during t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
colorbar
 
b2=b(:,:,2);
figure
%subplot(3,3,5)
pcolor(b2)
shading flat
h=colormap(jet);
rgb2hsv(h);
axis ij
%bmin=min(nonzeros(b1));
%bmax=max(max(b1));
%caxis([bmin bmax])
caxis([0.4 2])
title(['IBIs for 2nd and 3rd activations during t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
colorbar

figure
%g=subplot(3,3,7);
pcolor(c)
shading flat
colormap(jet);
rgb2hsv(h);
axis ij
caxis([-.3 .6])
title(['Change in IBIs for 2nd and 3rd beats during t1=',num2str(starttime),'s to t2=',num2str(endtime),'s']); 
colorbar
