function IBI=ibi(timeseries,thresh,samprate)
%this funciton computes interbeat intervals of a timeseries

lt=length(timeseries);
a=zeros(1,lt);
for i=2:lt-1;
    if timeseries(i)>thresh && timeseries(i-1)<timeseries(i)
       a(i)=i;
    end
end
b=nonzeros(a);   
c=b;
%f=0;
%minIBI=0.5*samprate;
%for j=2:length(b);
%    if (b(j)-b(j-1))<minIBI;
%        c(j-f)=[];
%        f=f+1
%    end
%end
IBI=diff(c);%/samprate;

%{
plot(IBI,'*')
ylabel('IBI (s)','FontSize',20)
xlabel('beat #','FontSize',20)
title('IBI distribution for filtered data','FontSize',20,'Interpreter','none');
%}