function [IBIdoc,timedoc]=IBIdoctor2(IBI,threshdoc,delay)
%this allows you to clump IBIs which are too short to the next IBI. it keeps clumping until
%the clump gets bigger then thresdoc.
%NOTE this assumes that the short IBIs after the real threshold crossing.
%this is almost always the case unless you set your upcross threshold too
%low.
%added delay (time of first crossing) to align with timeseries

IBIdoc1=zeros(1,length(IBI));

a=0;
for k=1:length(IBI)
    IBIdoc1(k)=IBI(k)+a;
        if IBIdoc1(k)<threshdoc
            a=IBIdoc1(k);
            IBIdoc1(k)=0;
        else
            a=0;
        end     
end

%truncate
IBIdoc=nonzeros(IBIdoc1);

%delay aligns with timeseries
timedoc=zeros(1,length(IBIdoc));
timedoc(1)=delay+IBIdoc(1);
for j=2:length(IBIdoc);
    timedoc(j)=IBIdoc(j)+timedoc(j-1);
end

