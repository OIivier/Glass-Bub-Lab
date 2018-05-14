function [IBIdoc,timedoc,meanIBIdoc,stdIBIdoc]=ibiBB(N,xcord,ycord,thresh,threshdoc)

timeseries=N(xcord,ycord,:);
lt=length(timeseries);

for i=1:length(timeseries);
     times(i)=timeseries(1,1,i);
     threshvec(i)=thresh;
    %to draw a line for thresh make threshvec constant for all times
end


sr=40; % 40Hz sampling
t=0:1/sr:lt/sr-1/sr; %lt length of time in secs


a=zeros(1,lt);

for i=1:lt-1;
    if (timeseries(i)<thresh && timeseries(i+1)>thresh)
       a(i)=(i/sr)+(thresh-timeseries(i))/(sr*(timeseries(i+1)-timeseries(i)));
    end
end
b=nonzeros(a);


IBI=diff(b);
meanIBI=mean(IBI);
stdIBI=std(IBI);

time(1)=IBI(1);
for j=2:length(IBI);
    time(j)=IBI(j)+sum(IBI(1:j-1));
end

[IBIdoc,timedoc]=IBIdoctor(IBI,threshdoc);

meanIBIdoc=mean(IBIdoc);
stdIBIdoc=std(IBIdoc);


