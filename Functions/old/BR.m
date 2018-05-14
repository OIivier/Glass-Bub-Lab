[v, ]=BR

%the Beleer-Reuter equations with diffusion. 
% c*dv/dt=*(I_s + I_Na +I_x1 + I_K1 + I_app) +D*(d^2/dx^2)v


c=

I_tot(2:lx-1,2:lx-1,1) = I_Na(2:lx1-1,2:lx2-1,1) + I_s(2:lx1-1,2:lx2-1,1) + I_x1(2:lx1-1,2:lx2-1,1) + I_K1(2:lx1-1,2:lx2-1,1) + I_app(2:lx1-1,2:lx2-1)
v(2:lx1-1,2:lx2-1,2) = v(2:lx1-1,2:lx2-1,1) + dt/c*(-I_tot(2:lx-1,2:lx-1) + (D/(dx^2))*(v(3:lx1,2:lx2-1,1)+v(1:lx1-2,2:lx2-1,1)+v(2:lx1-1,3:lx2,1)+v(2:lx1-1,1:lx2-2,1)-4*v(2:lx1-1,2:lx2-1,1)));

am(2:lx1-1,2:lx2-1,2)=am(2:lx-1,2:lx-1,1)-dt*(v(2:lx-1,2:lx-1,2)+47)./(exp(-0.1.*(v(2:lx-1,2:lx-1,2)+47))-1);
bm(2:lx-1,2:lx-1,2)=bm(2:lx-1,2:lx-1,1)+dt*(40*exp(-0.056*(v(2:lx-1,2:lx-1,2)+72)));
m(2:lx-1,2:lx-1,2)=m(2:lx-1,2:lx-1,1) +dt*(am(2:lx-1,2:lx-1,1).*(1-m(2:lx-1,2:lx-1,1))-bm(2:lx-1,2:lx-1,1).*m(2:lx-1,2:lx-1,1));
an(2:lx-1,2:lx-1,2)=an(2:lx-1,2:lx-1,1) +dt*(0.126*exp(-0.25*(v(2:lx-1,2:lx-1,2)+77)));
bn(2:lx-1,2:lx-1,2)=bn(2:lx-1,2:lx-1,1) +dt*1.7/(exp(-0.082*(v(2:lx-1,2:lx-1,2)+22.5))+1);
n(2:lx-1,2:lx-1,2)=n(2:lx-1,2:lx-1,1) +dt*(an(2:lx-1,2:lx-1,1).*(1-n(2:lx-1,2:lx-1,1))-bn(2:lx-1,2:lx-1,1).*n(2:lx-1,2:lx-1,1));
al(2:lx-1,2:lx-1,2)=al(2:lx-1,2:lx-1,1) +dt*(0.055*exp(-0.25*(v(2:lx-1,2:lx-1,2)+78))./(exp(-0.2*(v(2:lx-1,2:lx-1,2)+78))+1));
bl(2:lx-1,2:lx-1,2)=bl(2:lx-1,2:lx-1,1) +dt*(0.3./(exp(-0.1*(v(2:lx-1,2:lx-1,2)+32))+1));
l(2:lx-1,2:lx-1,2)=l(2:lx-1,2:lx-1,1) +dt*((al(2:lx-1,2:lx-1,1).*(1-l(2:lx-1,2:lx-1,1))-bl(2:lx-1,2:lx-1,1).*l(2:lx-1,2:lx-1,1)));

g_Na=23.0
g_Nac=0.003
Ena=?  %the eqm potential for Na
I_Na(2:lx-1,2:lx-1)=(gna*((m(2:lx-1,2:lx-1)).^3).*n(2:lx-1,2:lx-1).*l(2:lx-1,2:lx-1) + g_Nac).*(v(2:lx-1,2:lx-1)-Ena)