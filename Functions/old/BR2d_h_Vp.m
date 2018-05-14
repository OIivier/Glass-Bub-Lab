function [v,w]=FHN2d_h_Vp(n,time);

% solves 2D FHN eqns. with periodic boundary conditions using Euler method.
%FHN looks like:    dv/dt=-w-v*(v-a)*(v-1)+D*(d^2/dx^2)v+I
%                   dw/dt=e*(b*v-g*w-d)

% n is number of current sinks
% time is length of time

dx=0.02; % space step (cm), a myocyte is about 200 microns 
dt=0.1; % time step (ms)
% NB: D,dt,dx must be chosen such that 1-4*(dt*D/(dx^2))> = 0 to ensure
% proper averaging 
t=1:(time/dt)+1;
lt=length(t);

lx=201; %with dx-0.02cm this is about 4cmX4cm of tissue.
x1=1:lx; % this is the vertical coord.       
x2=1:lx;


ltI=100;  %how many time steps you wanna apply current for, here 10ms
I_app=zeros(lx,lx);  % applied current
I(:,2)=100; %in uA/cm^2
    
v=-70*ones(lx,lx,5);
c=(1e-3)*ones(lx,lx,2);
m=rand*ones(lx,lx,2);
n=rand*ones(lx,lx,2);
l=rand*ones(lx,lx,2);
d=rand*ones(lx,lx,2);
f=rand*ones(lx,lx,2);
z=rand*ones(lx,lx,2);    



cap=1;
D=0.001; % diffusion coeff. in units of cm^2/ms


p=randmat(n,lx);

for t=1:ltI; %the first loop defines equations with injected current
   
    I_tot(2:lx-1,2:lx-1,1) = I_Na(2:lx1-1,2:lx2-1,1) + I_s(2:lx1-1,2:lx2-1,1) + I_x1(2:lx1-1,2:lx2-1,1) + I_K1(2:lx1-1,2:lx2-1,1) - I_app(2:lx1-1,2:lx2)
    v(2:lx1-1,2:lx2-1,2) = v(2:lx1-1,2:lx2-1,1) + dt/cap*(-I_tot(2:lx-1,2:lx-1) + (D/(dx^2))*(v(3:lx1,2:lx2-1,1)+v(1:lx1-2,2:lx2-1,1)+v(2:lx1-1,3:lx2,1)+v(2:lx1-1,1:lx2-2,1)-4*v(2:lx1-1,2:lx2-1,1)));
     
    % because of periodic BC in the x1 coord, the space averaging must be explicitly defined at near the endpoints for v
    v(1,2:lx-1,2)=v(1,2:lx-1,1)+dt*(-w(1,2:lx-1,1)-v(1,2:lx-1,1).*(v(1,2:lx-1,1)-a).*(v(1,2:lx-1,1)-1)+(D/(dx^2))*(v(2,2:lx-1,1)+v(lx-1,2:lx-1,1)+v(1,3:lx,1)+v(1,1:lx-2,1)-4*v(1,2:lx-1,1))+I(1,2:lx-1));
    v(lx,:,2)=v(1,:,2);
    v(lx-1,2:lx-1,2)=v(lx-1,2:lx-1,1)+dt*(-w(lx-1,2:lx-1,1)-v(lx-1,2:lx-1,1).*(v(lx-1,2:lx-1,1)-a).*(v(lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(lx,2:lx-1,1)+v(lx-2,2:lx-1,1)+v(lx-1,3:lx,1)+v(lx-1,1:lx-2,1)-4*v(lx-1,2:lx-1,1))+I(lx-1,2:lx-1));
    % 1st order no-flux boundary conditions at ends of x2 
    v(:,1,2)=v(:,2,2);
    v(:,lx,2)=v(:,lx-1,2);
    
    
am(2:lx1-1,2:lx2-1,2)=-(v(2:lx-1,2:lx-1,2)+47)./(exp(-0.1.*(v(2:lx-1,2:lx-1,2)+47))-1);
bm(2:lx-1,2:lx-1,2)= (40*exp(-0.056*(v(2:lx-1,2:lx-1,2)+72)));
m(2:lx-1,2:lx-1,2)=m(2:lx-1,2:lx-1,1) +dt*(am(2:lx-1,2:lx-1,1).*(1-m(2:lx-1,2:lx-1,1))-bm(2:lx-1,2:lx-1,1).*m(2:lx-1,2:lx-1,1));
an(2:lx-1,2:lx-1,2)=(0.126*exp(-0.25*(v(2:lx-1,2:lx-1,2)+77)));
bn(2:lx-1,2:lx-1,2)=1.7/(exp(-0.082*(v(2:lx-1,2:lx-1,2)+22.5))+1);
n(2:lx-1,2:lx-1,2)=n(2:lx-1,2:lx-1,1) +dt*(an(2:lx-1,2:lx-1,1).*(1-n(2:lx-1,2:lx-1,1))-bn(2:lx-1,2:lx-1,1).*n(2:lx-1,2:lx-1,1));
al(2:lx-1,2:lx-1,2)=(0.055*exp(-0.25*(v(2:lx-1,2:lx-1,2)+78))./(exp(-0.2*(v(2:lx-1,2:lx-1,2)+78))+1));
bl(2:lx-1,2:lx-1,2)=(0.3./(exp(-0.1*(v(2:lx-1,2:lx-1,2)+32))+1));
l(2:lx-1,2:lx-1,2)=l(2:lx-1,2:lx-1,1) +dt*((al(2:lx-1,2:lx-1,1).*(1-l(2:lx-1,2:lx-1,1))-bl(2:lx-1,2:lx-1,1).*l(2:lx-1,2:lx-1,1)));

g_Na=4.0;
g_Nac=0.003;
Ena=50.0  %the eqm potential for Na
I_Na(2:lx-1,2:lx-1,2)=I_Na(2:lx-1,2:lx-1,1) +dt*(g_Na*((m(2:lx-1,2:lx-1)).^3).*n(2:lx-1,2:lx-1).*l(2:lx-1,2:lx-1) + g_Nac).*(v(2:lx-1,2:lx-1)-Ena)


af(2:lx-1,2:lx-1,2)= 0.012*exp(0.008*(v(2:lx-1,2:lx-1,1)+28))./(exp(0.15*(v(2:lx-1,2:lx-1,1)+28))+1);
bf(2:lx-1,2:lx-1,2)= 0.0065*exp(-0.02*(v(2:lx-1,2:lx-1,1)+30))./(exp(-0.2*(v(2:lx-1,2:lx-1,1)+30))+1);
f(2:lx-1,2:lx-1,2)= f(2:lx-1,2:lx-1,1) +dt*((af(2:lx-1,2:lx-1,1).*(1-f(2:lx-1,2:lx-1,1))-(bf(2:lx-1,2:lx-1,1).*f(2:lx-1,2:lx-1,1))));
ad(2:lx-1,2:lx-1,2)= 0.095*exp(-0.01*(v(2:lx-1,2:lx-1,1)-5))./(exp(-0.072(v(2:lx-1,2:lx-1,1)-5))+1);
bd(2:lx-1,2:lx-1,2)= 0.07*exp(-0.017*(v(2:lx-1,2:lx-1,1)+44))./(exp(0.05*(v(2:lx-1,2:lx-1,1)+44))+1);
d(2:lx-1,2:lx-1,2)= d(2:lx-1,2:lx-1,1) +dt*((ad(2:lx-1,2:lx-1,1).*(1-d(2:lx-1,2:lx-1,1))-(bd(2:lx-1,2:lx-1,1).*d(2:lx-1,2:lx-1,1))));

g_Ca=0.09;
c(2:lx-1,2:lx-1,2)=c(2:lx-1,2:lx-1,1) +dt*(0.07*((1e-4)-c(2:lx-1,2:lx-1,1))-(10^-4)*I_na(2:lx-1,2:lx-1,1));
Cai(2:lx-1,2:lx-1,1)=(1e-3)*c(2:lx-1,2:lx-1,1);
I_Ca(2:lx-1,2:lx-1,2)=I_Ca(2:lx-1,2:lx-1,1) +dt*(g_Ca.*f(2:lx-1,2:lx-1,1).*d(2:lx-1,2:lx-1,1).*(v(2:lx-1,2:lx-1,1)+82.3+13.0287*log(Cai(2:lx-1,2:lx-1,1))))


I_K(2:lx-1,2:lx-1,2)=1.4*(exp(0.04*(v(2:lx-1,2:lx-1,1)+85))-1)./(exp(0.08*(v(2:lx-1,2:lx-1,1)+53))+exp(0.04*(v(2:lx-1,2:lx-1,1)+53)));

az(2:lx-1,2:lx-1,1)=0.0005*exp(0.083*(v(2:lx-1,2:lx-1,1)+50))./(exp(0.057*(v(2:lx-1,2:lx-1,1)+50))+1);
bz(2:lx-1,2:lx-1,1)=0.0013*(exp(-0.06*(v(2:lx-1,2:lx-1,1)+20)))./(exp(-0.04*(v(2:lx-1,2:lx-1,1)+20))+1);
z(2:lx-1,2:lx-1,2)=z(2:lx-1,2:lx-1,1) +dt*((az(2:lx-1,2:lx-1,1).*(1-z(2:lx-1,2:lx-1,1))-(bz(2:lx-1,2:lx-1,1).*z(2:lx-1,2:lx-1,1))));
I_z(2:lx-1,2:lx-1,2)=z(2:lx-1,2:lx-1,1)*0.8*(exp(0.04*(v(2:lx-1,2:lx-1,2)+77))-1)./exp(0.04*(v(2:lx-1,2:lx-1,1)+35))
    
    % BC (periodic at ends of x2, no-flux at ends of x1) for w

    
    for i1=1:lx; %every time we update v, w the current sink positions must be rezeroed 
        for j1=1:lx;
        
          if p(i1,j1)==0
            v(i1,j1,2)=p(i1,j1);
          end
        end 
    end
    
    v(:,:,3)=v(:,:,2); % want to save the end of pulse for viewing
    
    
    v(:,:,1)=v(:,:,2); %this shifts updated time matriz to previous time matriz spot for next iteration
    m(:,:,1)=m(:,:,2);
    n(:,:,1)=n(:,:,2);
    l(:,:,1)=l(:,:,2);
    f(:,:,1)=f(:,:,2);
    d(:,:,1)=d(:,:,2);
    c(:,:,1)=c(:,:,2);
    z(:,:,1)=z(:,:,2);
    
    
end


for t=ltI:time/dt; %2nd for loop is identical to 1st except without injected current. also here we calcualte Vp
    
    

    v(2:lx-1,2:lx-1,2)=v(2:lx-1,2:lx-1,1)+dt*(-w(2:lx-1,2:lx-1,1)-v(2:lx-1,2:lx-1,1).*(v(2:lx-1,2:lx-1,1)-a).*(v(2:lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(3:lx,2:lx-1,1)+v(1:lx-2,2:lx-1,1)+v(2:lx-1,3:lx,1)+v(2:lx-1,1:lx-2,1)-4*v(2:lx-1,2:lx-1,1)));
     
    % because of periodic BC in the x1 coord, the space averaging must be explicitly defined at near the endpoints for v
    v(1,2:lx-1,2)=v(1,2:lx-1,1)+dt*(-w(1,2:lx-1,1)-v(1,2:lx-1,1).*(v(1,2:lx-1,1)-a).*(v(1,2:lx-1,1)-1)+(D/(dx^2))*(v(2,2:lx-1,1)+v(lx-1,2:lx-1,1)+v(1,3:lx,1)+v(1,1:lx-2,1)-4*v(1,2:lx-1,1)));
    v(lx,:,2)=v(1,:,2);
    v(lx-1,2:lx-1,2)=v(lx-1,2:lx-1,1)+dt*(-w(lx-1,2:lx-1,1)-v(lx-1,2:lx-1,1).*(v(lx-1,2:lx-1,1)-a).*(v(lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(lx,2:lx-1,1)+v(lx-2,2:lx-1,1)+v(lx-1,3:lx,1)+v(lx-1,1:lx-2,1)-4*v(lx-1,2:lx-1,1)));
    % 1st order no-flux boundary conditions at ends of x2 NB to get 2nd
    % order approx (less error) must generate and solve for phantom points
    % -see notes
    v(:,1,2)=v(:,2,2);
    v(:,lx,2)=v(:,lx-1,2);
    
    w(1:lx-1,:,2)=w(1:lx-1,:,1)+dt*e*(b*v(1:lx-1,:,1)-g*w(1:lx-1,:,1)-d);
    
    % BC (periodic at ends of x2, no-flux at ends of x1) for w
    w(lx,:,2)=w(1,:,2); 
    w(:,1,2)=w(:,2,2);
    w(:,lx,2)=w(:,lx-1,2);
    
    
    %more snapshots to store
    if t==999;    
        v(:,:,4)=v(:,:,2);
        w(:,:,4)=w(:,:,2);
    elseif t==1999;
        v(:,:,5)=v(:,:,2);
        w(:,:,5)=w(:,:,2);
    end

    for i1=1:lx;
        for j1=1:lx;
        
          if p(i1,j1)==0
            v(i1,j1,2)=p(i1,j1);
            w(i1,j1,2)=p(i1,j1);
          end
        end 
    end
    
    
    %this is the propagation velocity, Vp, estimator
    if v(20,50,2)>=0.3 && v(20,50,1)<0.3 %we catch v on the upstroke and make sure it's not a sink.
       t1(t)=t;%this sets up so that we keep only first time that this happens.
    else t1(t)=0;
    end 
    if v(20,200,2)>=0.3 && v(20,200,1)<0.3
       t2(t)=t;
    else t2(t)=0;
    end
    
    
    v(:,:,1)=v(:,:,2); %this shifts updated time matriz to previous time matriz spot for next iteration
    w(:,:,1)=w(:,:,2);
end %end of time loop
    
if isempty(find(t1,1, 'first'))==0 &&  isempty(find(t2,1, 'first'))==0
    t11=find(t1,1, 'first')
    t22=find(t2,1, 'first')
    Vp=(150*0.02)/(0.1*(t22-t11))*1e3 %prop vel. in cm/s
else Vp=0
end


newplot
subplot(2,3,1); pcolor(v(:,:,3)); %plot(1:lx,v(10,:,3)); %; %check out in space @ timesteps ______
caxis([-1 1])
shading interp
colormap(jet)
subplot(2,3,2); pcolor(v(:,:,4)); %plot(1:lx,v(10,:,4)); %
caxis([-1 1])
shading interp
colormap(jet)
subplot(2,3,3); pcolor(v(:,:,5)); %plot(1:lx,v(10,:,5)); %;
caxis([-1 1])
shading interp
colormap(jet)
subplot(2,3,4); pcolor(w(:,:,3)); %plot(1:lx,w(10,:,3)); %
caxis([-1 1])
shading interp
colormap(jet)
subplot(2,3,5); pcolor(w(:,:,4)); %plot(1:lx,w(10,:,4)); %
caxis([-1 1])
shading interp
colormap(jet)
subplot(2,3,6); pcolor(w(:,:,5)); %plot(1:lx,w(10,:,5)); %
caxis([-1 1])
shading interp
colormap(jet)
colorbar('north')
