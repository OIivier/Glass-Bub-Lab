function [p,time,ibi_s,ibi_f]=phase_chng_simple(s1,s2,s3,dt)
%this measures the plahse change of the slow oscillator, in absence of
%another oscillator (s1), and in its presence (s2).
%s3 is the timeseries of the fast pacemaker site.

ls=length(s1);
% length of s1, and s2 should be same

a=zeros(1,ls);
b=zeros(1,ls);
c=zeros(1,ls);

for i=2:ls-1;
    if s1(i)<s1(i-1) && s1(i)<s1(i+1)
       a(i)=i;
    end
end

%make sure all ibis of s1 are equal when there is no fast pacemaker
ibi_s=diff(nonzeros(a))*dt;
%if that's true just use one of them as representing the unpertrubed period
ibi_s1=ibi_s(10);

%now we see check the IBI of the slow pacemaker in the presence of the fast pacemaker
for i=2:ls-1;
    if s2(i)<s2(i-1) && s2(i)<s2(i+1)
       b(i)=i;
    end
end

b1=nonzeros(b);
lb1=length(b1);
time=zeros(1,lb1-1);
p=zeros(1,lb1-1);

for i=1:lb1-1
    p(i)=(b1(i+1)-(b1(i)+ibi_s1))/ibi_s1;
    time(i)=b1(i+1)*dt;
end

%find the preiod of the fast oscillator
for i=2:ls-1;
    if s3(i)<s3(i-1) && s3(i)<s3(i+1)
       c(i)=i;
    end
end

%check all ibis of s3 are equal
ibi_f=diff(nonzeros(c))*dt;
%if that's true just use one of them as representing the unpertrubed period
ibi_f1=ibi_f(10);

tot_time=0:dt:ls*dt-dt;

figure
plot(tot_time,s1)
hold
plot(tot_time,s2,'g')
plot(tot_time,s3,'r')
xlabel('time[s]')
ylabel('v(t)')
title(['time series of 3 pacemaker sites for IBI_s=',num2str(ibi_s1),'s, IBI_f=', num2str(ibi_f1),'s'])  
legend('SP with no FP','SP with FP','FP')
hold

figure
plot(time,p,'*-')
xlabel('time[s]')
ylabel('phase shift [%]')
title(['Phase shifting of slow pacemaker for IBI_s=',num2str(ibi_s1),'s, IBI_f=', num2str(ibi_f1),'s'])  