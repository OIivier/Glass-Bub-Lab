%plot Vp vs rho (density iof sinks) for different D

clear
tic
spaces=15;
n=[0,100,200,linspace(250,700,10)]; %vector of heterogeneities numbers
ltime=4000;

%vp1=zeros(spaces,1);
%vp2=zeros(spaces,1);
vp6=zeros(spaces,1);
%vp4=zeros(spaces,1);
%vp4=zeros(spaces,1);

D=[0.00001,0.00035,0.0007,0.001, 0.0014];


for i=1:length(n);
%    [vp1(i)]=FHN2d_sinks(ltime,D(1),n(i));
%    [vp2(i)]=FHN2d_sinks(ltime,D(2),n(i));
     [vp6(i)]=FHN2d_sinks(ltime,D(4),0.02,0.01, n(i));
%    [vp4(i)]=FHN2d_sinks(ltime,D(4),n(i));
%    [vp3(5)]=FHN2d_sinks(ltime,D(5),n(i));
end

%subplot(2,1,1);errorbar(rho,Vp4,Vp4err,'.')
%subplot(2,1,2);plot(rho,lVp4,'*')
%subplot(2,1,3);plot(rho,Vp4err,'--')
%hist(Vp4)
%hold on
%plot(n/(200*0.02)^2,Vp2,'b')
%hold on
%plot(n/(200*0.02)^2,Vp3,'g')
%hold on
%plot(n/(200*0.02)^2,Vp4,'m')
%hold on
%plot(n/(200*0.02)^2,Vp5,'r')



%title('Effect of Current Sink Density on Effect of Current Sink Density on Mean Wave Propagation Velocity')
%xlabel('sink density [n/cm^2]')
%ylabel('propagation velocity [cm/s]')
%legend('D=0.0005cm^2/ms','D=0.0007cm^2/ms', 'D=0.0011cm^2/ms')

%h = legend('D=0.0003 cm^2/s','D=0.0005 cm^2/s','D=0.0007 cm^2/s', 'D=0.0010 cm^2/s', 'D=0.0012 cm^2/s', 5);
%set(h,'Interpreter','none')

toc
save('c_vs_n_a02_e01_D01')
