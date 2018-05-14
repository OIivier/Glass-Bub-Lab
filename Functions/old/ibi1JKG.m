function [IBI,time]=ibi1JKG(N,xcord,ycord,thresh)
%this funciton computes interbeat intervals of a timeseries
%Make sure you truncate your control
close all;


sr=40;
timeseries=N(xcord,ycord,:);
lt=length(timeseries);

for i=1:length(timeseries);
     times(i)=timeseries(1,1,i);
     threshvec(i)=thresh;
end


t=0:1/sr:lt/sr-1/sr;
figure(1);
subplot(2,1,1);
hold on;
plot(t,times,'b');
plot(t,threshvec,'r');
xlabel('Time (s)','Fontsize',16);
ylabel('Light intensity','Fontsize',16);
title('Time Series','Fontsize',20);
legend('Timeseries','Threshold','Fontsize',12);

a=zeros(1,lt);
sr=40; % This is always the case

for i=1:lt-1;
    if (timeseries(i)<thresh && timeseries(i+1)>thresh)
       a(i)=(i/sr)+(thresh-timeseries(i))/(sr*(timeseries(i+1)-timeseries(i)));
    end
end
b=nonzeros(a);


IBI=diff(b);

time(1)=IBI(1);
for j=2:length(IBI);
    time(j)=IBI(j)+sum(IBI(1:j-1));
end

% IBIdoc(1)=IBI(1);
% for k=2:length(IBI);
%     if IBI(k)<thresh2;
%         IBIdoc(k)=IBI(k)+IBI(k-1);
%     end
%     if IBI(k)>thresh2
%         IBIdoc(k)=IBI(k);
%     end
%     
% end

% subplot(3,1,3);
% plot(IBI,'r*');
% ylim([0 2]);
% ylabel('IBI (s)','FontSize',20)
% xlabel('beat #','FontSize',20)
% title('IBI distribution for filtered data','FontSize',20,'Interpreter','none');
% ylim([0.5 1.5]);

subplot(2,1,2);
hold on;
plot(time,IBI,'r*');  title('Time aligned IBIs','FontSize',20); ylim([0.5 2]);
ylabel('IBI (s)','FontSize',16)
xlabel('Beat Number','FontSize',16);
 
% % subplot(3,1,3);
% plot(time,IBIdoc,'g*'); 
% title('Doctored IBIs','FontSize',20); ylim([0.5 1.5]);
% ylabel('IBI (s)','FontSize',16)
% xlabel('Beat Number','FontSize',16)




