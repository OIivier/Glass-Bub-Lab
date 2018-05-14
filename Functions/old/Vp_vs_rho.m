n=0:250:6000;
vpm=zeros(25,1);
s=zeros(25,1);
l=zeros(25,1);

for i=1:25
    [vpm(i),n(i)]=FHN2d_breaks100(n(i));
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
