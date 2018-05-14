function [IBI,IBIdoc,time,timedoc,meanIBI,meanIBIdoc,stdIBI,stdIBIdoc]=IBIplot_manythresh(N,xcord,ycord,threshs,threshdoc,switchtimes)

%this funciton computes interbeat intervals of 9 pixels and displays their
%statistics depending of how you set the threshold for getting rid of short
%IBIs (doctoring the IBIs)

%original IBI detector made by Bart in 2008
% Jake Gabriels coded the parts that align times of IBIs march 2009
%Bart added the statistical plots, lt, threshdoc april 2009

timeseries=N(ycord,xcord,:);
 
lt=length(timeseries);
for i=1:length(timeseries);
     times(i)=timeseries(1,1,i);
end

% To make vertical lines on the graphs whenever there is a change in the dynamics
vertime = (min(timeseries):2:max(timeseries));
vertIBI = 0:0.1:2.5;

sr=40; % 40Hz sampling
t=0:1/sr:lt/sr-1/sr; %lt length of time in secs

if numel(switchtimes)==0
    switchtimes=lt/40;
else
    switchtimes=[switchtimes lt/40];
    threshs = [threshs threshs(1)];
end


newt(1:switchtimes(1)*40,1)=t(1:switchtimes(1)*40)';
threshvec(1:switchtimes(1)*40,1)=threshs(1)';

if length(switchtimes)~=1
    for j=2:length(threshs)
        newt(switchtimes(j-1)*40:switchtimes(j)*40,j)=t(switchtimes(j-1)*40:switchtimes(j)*40)';
        threshvec(switchtimes(j-1)*40:switchtimes(j)*40,j)=threshs(j)';
    end
end

figure
subplot(2,1,1);
hold on;
plot(t,times,'b');
plot(newt(1:switchtimes(1)*40,1),threshvec(1:switchtimes(1)*40,1),'r');

if length(switchtimes)~=1
    for j=2:length(threshs);
        plot(newt(switchtimes(j-1)*40:switchtimes(j)*40,j),threshvec(switchtimes(j-1)*40:switchtimes(j)*40,j),'r');
    end
end

for i=1:length(switchtimes);   
    plot(switchtimes(i),vertime,'c.','MarkerSize',5);
end

xlabel('Time (s)','Fontsize',10);
ylabel('Zero-meaned Light intensity','Fontsize',10);
title(['Time Series for Pixel (',num2str(xcord),',',num2str(ycord),')'],'FontSize',12,'Interpreter','none');
%legend('Timeseries','Threshold','Fontsize',8);

a=zeros(lt,length(switchtimes));

for i=1:switchtimes(1)*40-1;
    if (timeseries(i)<threshs(1) && timeseries(i+1)>threshs(1))
        a(i)=(i/sr)+(threshs(1)-timeseries(i))/(sr*(timeseries(i+1)-timeseries(i)));
    end
end

if length(switchtimes)~=1
    for j=2:length(threshs);
        for i=switchtimes(j-1)*40:switchtimes(j)*40-1;
            if (timeseries(i)<threshs(j) && timeseries(i+1)>threshs(j))
                a(i)=(i/sr)+(threshs(j)-timeseries(i))/(sr*(timeseries(i+1)-timeseries(i)));
            end
        end
    end
end

b=nonzeros(a);

IBI=diff(b);
meanIBI=mean(IBI);
stdIBI=std(IBI);

%this will always be a little delayed relative to t since the first IBI
%doesn't happen at t=0..added b(1) april 29 call ibiBBplot1.m

%error if no crossing times... numel(b)=0
%thus

if numel(b)==0
    b=0;
end

time(1)=b(1)+IBI(1);
for j=2:length(IBI);
    time(j)=IBI(j)+time(j-1);
end

[IBIdoc,timedoc]=IBIdoctor2(IBI,threshdoc,b(1));

meanIBIdoc=mean(IBIdoc);
stdIBIdoc=std(IBIdoc);

subplot(2,1,2);
hold on;
plot(time,IBI,'r*'); 
plot(timedoc,IBIdoc,'g*');
for i=1:length(switchtimes)
    plot(switchtimes(i),vertIBI,'c.','MarkerSize',5);
end
title('Interbeat Intervals','FontSize',12); 
ylim([0 2.5]);
ylabel('IBI (s)','FontSize',10);
xlabel('time (s)','FontSize',10);
 



