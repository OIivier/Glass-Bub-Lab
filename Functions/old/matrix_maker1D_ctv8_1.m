function [f1,f2]=matrix_maker1D_ctv8_1(f);

lx=6; % how many sequential IBIs the code outputs for each stimulus phase.

f=f(:,2); %take out the crossing time vector

f1=zeros(lx,length(f)/(lx));

for i=1:length(f)/lx;
    f1(:,i)=f(lx*(i-1)+1:lx*i);   
end

%normalize the first 3 sets of affected ibis, given that f(2,1) is
%an unperturbed ibi.

f2=f1;
f2(3,:)=f1(3,:)/f1(2,2);
f2(4,:)=(f1(3,:)+f1(4,:))/f1(2,2);
f2(5,:)=(f1(5,:)+f1(4,:)+f1(3,:))/f1(2,2);

%assuming that your stimulus covers the entire cycle exactly, make a vector
%of phases,
phase=linspace(0,1,length(f1(1,:)));

figure
plot(phase,f2(3,:),'*')
hold
plot(phase,f2(4,:),'g*')
plot(phase,f2(5,:),'r*')
ylim([0,3.5]);
