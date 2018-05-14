function [phase,f1,f2]=matrix_maker1D_ctIBI(f,phases);

lx=6; % how many sequential IBIs the code outputs for each stimulus phase.

f=f(:,3); %take out the crossing time vector

f1=zeros(lx,length(f)/(lx));

for i=1:length(f)/lx;
    f1(:,i)=f(lx*(i-1)+1:lx*i);   
end

%normalize the first 3 sets of affected ibis, given that f(2,1) is
%an unperturbed ibi.

f2=f1;
f2(4,:)=f1(4,:)/f1(3,2);
f2(5,:)=(f1(5,:)+f1(4,:))/f1(3,2);
f2(6,:)=(f1(6,:)+f1(5,:)+f1(4,:))/f1(3,2);

%assuming that your stimulus covers the entire cycle exactly, make a vector
%of phases,T

%phase=linspace(0.30,0.36,length(f1(1,:)));
phase=phases(:,3);

figure
plot(phase,f2(4,:),'*')
hold
plot(phase,f2(5,:),'g*')
plot(phase,f2(6,:),'r*')
ylim([0,3.5]);
xlim([0,1]);	
xlabel(texlabel('phi'))
ylabel('T_i/T_0')
title('phase resetting at a distance of Xcm')
