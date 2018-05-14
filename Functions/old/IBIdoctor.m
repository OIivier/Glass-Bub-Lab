function [IBIdoc,timedoc]=IBIdoctor(IBI,threshdoc)
%this allows you to clump IBIs which are too short
%NOTE: two consecutive short IBIs will not be grouped
IBIdoc=zeros(length(IBI));
IBIdoc(1)=IBI(1);

l=0;
m=0;
for k=1:length(IBI)-1
     if (IBI(k)<threshdoc)
        IBIdoc(k-m)=IBI(k)+IBI(k+1);
        l=l+1;
        m=m+1;
     else 
         IBIdoc(k-m+l)=IBI(k);
         l=0;
     end
end

%truncate
IBIdoc=nonzeros(IBIdoc);

%assumes you start at first threshold crossing of timeseries.
timedoc(1)=IBIdoc(1);
for j=2:length(IBIdoc);
    timedoc(j)=IBIdoc(j)+sum(IBIdoc(1:j-1));
end

