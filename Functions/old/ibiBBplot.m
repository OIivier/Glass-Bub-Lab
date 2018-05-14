function [IBI,IBIdoc,time,timedoc,meanIBI,meanIBIdoc,stdIBI,stdIBIdoc]=ibiBBplot(N,xcord,ycord,thresh,threshdoc)
%this funciton computes interbeat intervals of 9 pixels and displays their
%statistics depending of how you set the threshold for getting rid of short
%IBIs (doctoring the IBIs)

%original IBI detector made by Bart in 2008
% Jake Gabriels coded the parts that align times of IBIs march 2009
%Bart added the statistical plots, lt, threshdoc april 2009

timeseries=N(xcord,ycord,:);
lt=length(timeseries);

for i=1:length(timeseries);
     times(i)=timeseries(1,1,i);
     threshvec(i)=thresh;
    %to draw a line for thresh make threshvec constant for all times
end


sr=40; % 40Hz sampling
t=0:1/sr:lt/sr-1/sr; %lt length of time in secs

figure

subplot(2,1,1);
hold on;
plot(t,times,'b');
plot(t,threshvec,'r');
xlabel('Time (s)','Fontsize',10);
ylabel('Light intensity','Fontsize',10);
title(['Time Series for Pixel (',num2str(xcord),',',num2str(ycord),')'],'FontSize',12,'Interpreter','none');
%legend('Timeseries','Threshold','Fontsize',8);

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

%this will always be a little delayed relative to t since the first IBI
%doesn't happen at t=0..
time(1)=IBI(1);
for j=2:length(IBI);
    time(j)=IBI(j)+sum(IBI(1:j-1));
end

[IBIdoc,timedoc]=IBIdoctor(IBI,threshdoc);

meanIBIdoc=mean(IBIdoc);
stdIBIdoc=std(IBIdoc);

subplot(2,1,2);
hold on;
plot(time,IBI,'r*'); 
plot(timedoc,IBIdoc,'g*');
title('Interbeat Intervals','FontSize',12); 
ylim([0 2.5]);
ylabel('IBI (s)','FontSize',10);
xlabel('time (s)','FontSize',10);
 



