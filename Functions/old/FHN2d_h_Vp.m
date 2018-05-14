function [v,w,Vp]=FHN2d_h_Vp(n,time);

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
I=zeros(lx,lx);  % applied current
    I(:,2)=1; 
    
v=zeros(lx,lx,5); %these also define initial conditions
w=zeros(lx,lx,5);



D=0.0007; % diffusion coeff. in units of cm^2/ms
a=0.02;
e=0.01;
b=0.5;
g=1;
d=0;

p=randmat(n,lx);

for t=1:ltI; %the first loop defines equations with injected current
   
    v(2:lx-1,2:lx-1,2)=v(2:lx-1,2:lx-1,1)+dt*(-w(2:lx-1,2:lx-1,1)-v(2:lx-1,2:lx-1,1).*(v(2:lx-1,2:lx-1,1)-a).*(v(2:lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(3:lx,2:lx-1,1)+v(1:lx-2,2:lx-1,1)+v(2:lx-1,3:lx,1)+v(2:lx-1,1:lx-2,1)-4*v(2:lx-1,2:lx-1,1))+I(2:lx-1,2:lx-1));
     
    % because of periodic BC in the x1 coord, the space averaging must be explicitly defined at near the endpoints for v
    v(1,2:lx-1,2)=v(1,2:lx-1,1)+dt*(-w(1,2:lx-1,1)-v(1,2:lx-1,1).*(v(1,2:lx-1,1)-a).*(v(1,2:lx-1,1)-1)+(D/(dx^2))*(v(2,2:lx-1,1)+v(lx-1,2:lx-1,1)+v(1,3:lx,1)+v(1,1:lx-2,1)-4*v(1,2:lx-1,1))+I(1,2:lx-1));
    v(lx,:,2)=v(1,:,2);
    v(lx-1,2:lx-1,2)=v(lx-1,2:lx-1,1)+dt*(-w(lx-1,2:lx-1,1)-v(lx-1,2:lx-1,1).*(v(lx-1,2:lx-1,1)-a).*(v(lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(lx,2:lx-1,1)+v(lx-2,2:lx-1,1)+v(lx-1,3:lx,1)+v(lx-1,1:lx-2,1)-4*v(lx-1,2:lx-1,1))+I(lx-1,2:lx-1));
    % 1st order no-flux boundary conditions at ends of x2 
    v(:,1,2)=v(:,2,2);
    v(:,lx,2)=v(:,lx-1,2);
    
    w(1:lx-1,:,2)=w(1:lx-1,:,1)+dt*e*(b*v(1:lx-1,:,1)-g*w(1:lx-1,:,1)-d);
    
    % BC (periodic at ends of x2, no-flux at ends of x1) for w
    w(lx,:,2)=w(1,:,2); 
    w(:,1,2)=w(:,2,2);
    w(:,lx,2)=w(:,lx-1,2);
    
    for i1=1:lx; %every time we update v, w the current sink positions must be rezeroed 
        for j1=1:lx;
        
          if p(i1,j1)==0
            v(i1,j1,2)=p(i1,j1);
            w(i1,j1,2)=p(i1,j1);
          end
        end 
    end
    
    v(:,:,3)=v(:,:,2); % want to save the end of pulse for viewing
    v(:,:,3)=v(:,:,2);
    
    v(:,:,1)=v(:,:,2); %this shifv(5,15,2)ts updated time matriz to previous time matriz spot for next iteration
    w(:,:,1)=w(:,:,2);
    
   
    
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
