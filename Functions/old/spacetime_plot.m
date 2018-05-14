function f1=spacetime_plot(spacetime)

lx=20

f1=matrix_maker1D(spacetime,lx);
t=linspace(5500+0.05,6000,length(f1(:,1)));
x=linspace(1,lx,length(f1(1,:)));
surf(x,t,f1)
shading interp
view(0,90)
axis ij
xlim([min(x), max(x)]);
ylim([min(t), max(t)]);
colorbar
xlabel('cell #')
ylabel('time [ms]')
title(['Spacetime plot of pacing from distance of 0.4cm at 0.3*T'], 'Fontsize',12)
%title(['Spacetime plot of resetting failure from distance of 0.4cm at ',texlabel('phi'),'=0.1667'], 'Fontsize',12)