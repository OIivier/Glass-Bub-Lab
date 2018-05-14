function [v,w]=FHN2d_h(rho,D,lx1,time);

% solves 2D FHN eqns. with periodic boundary conditions using Euler method.
%FHN looks like:    dv/dt=-w-v*(v-a)*(v-1)+D*(d^2/dx^2)v+I
%                   dw/dt=e*(b*v-g*w-d)

% rho is density of current sinks [#/cm^2]
% time is length of time (ms)
%lx1 is number of rows of space steps (remember to include the extra space
%for periodic BC)

dt=0.1; % time step (ms)
% NB: D,dt,dx must be chosen such that 1-4*(dt*D/(dx^2))> = 0 to ensure
% proper averaging 
t=1:(time/dt)+1;
lt=length(t);
dx=0.02; % space step (cm), a myocyte is about 200 microns 
 % this is the horizontal coord.     
lx2=201; %the extra one is because the no-flux BC



ltI=100;  %how many time steps you wanna apply current for, here 10ms
I=zeros(lx1,lx2);  % applied current
    I(:,2)=1; 
    
v=zeros(lx1,lx2,5); %these also define initial conditions
w=zeros(lx1,lx2,5);



%D=0.0007; % diffusion coeff. in units of cm^2/ms
a=0.02;
e=0.01;
b=0.5;
g=1;
d=0;

n=rho*(lx1-1)*lx2*(dx^2);
p=randmat2(n,lx1-1,lx2);

for t=1:ltI; %the first loop defines equations with injected current
      
    v(2:lx1-1,2:lx2-1,2)=v(2:lx1-1,2:lx2-1,1)+dt*(-w(2:lx1-1,2:lx2-1,1)-v(2:lx1-1,2:lx2-1,1).*(v(2:lx1-1,2:lx2-1,1)-a).*(v(2:lx1-1,2:lx2-1,1)-1)+(D/(dx^2))*(v(3:lx1,2:lx2-1,1)+v(1:lx1-2,2:lx2-1,1)+v(2:lx1-1,3:lx2,1)+v(2:lx1-1,1:lx2-2,1)-4*v(2:lx1-1,2:lx2-1,1))+I(2:lx1-1,2:lx2-1));
     
    % because of periodic BC in the x1 coord, the space averaging must be explicitly defined at near the endpoints for v
    v(1,2:lx2-1,2)=v(1,2:lx2-1,1)+dt*(-w(1,2:lx2-1,1)-v(1,2:lx2-1,1).*(v(1,2:lx2-1,1)-a).*(v(1,2:lx2-1,1)-1)+(D/(dx^2))*(v(2,2:lx2-1,1)+v(lx1-1,2:lx2-1,1)+v(1,3:lx2,1)+v(1,1:lx2-2,1)-4*v(1,2:lx2-1,1))+I(1,2:lx2-1));
    v(lx1,:,2)=v(1,:,2);
    v(lx1-1,2:lx2-1,2)=v(lx1-1,2:lx2-1,1)+dt*(-w(lx1-1,2:lx2-1,1)-v(lx1-1,2:lx2-1,1).*(v(lx1-1,2:lx2-1,1)-a).*(v(lx1-1,2:lx2-1,1)-1)+(D/(dx^2))*(v(lx1,2:lx2-1,1)+v(lx1-2,2:lx2-1,1)+v(lx1-1,3:lx2,1)+v(lx1-1,1:lx2-2,1)-4*v(lx1-1,2:lx2-1,1))+I(lx1-1,2:lx2-1));
    % 1st order no-flux boundary conditions at ends of x2 
    v(:,1,2)=v(:,2,2);
    v(:,lx2,2)=v(:,lx2-1,2);
    
    w(1:lx1-1,:,2)=w(1:lx1-1,:,1)+dt*e*(b*v(1:lx1-1,:,1)-g*w(1:lx1-1,:,1)-d);
    
    % BC (periodic at ends of x2, no-flux at ends of x1) for w
    w(lx1,:,2)=w(1,:,2); 
    w(:,1,2)=w(:,2,2);
    w(:,lx2,2)=w(:,lx2-1,2);
    
    for x1=1:lx1-1; %every time we update v, w the current sink positions must be rezeroed 
        for x2=1:lx2;
          if p(x1,x2)==0
            v(x1,x2,2)=0;
            w(x1,x2,2)=0;
          end
        end 
    end
    
    v(:,:,3)=v(:,:,2); % want to save the end of pulse for viewing
    v(:,:,3)=v(:,:,2);
    
    v(:,:,1)=v(:,:,2); %this shifv(5,15,2)ts updated time matriz to previous time matriz spot for next iteration
    w(:,:,1)=w(:,:,2);
    
end


for t=ltI+1:time/dt; %2nd for loop is identical to 1st except without injected current. also here we calcualte Vp
    
    v(2:lx1-1,2:lx2-1,2)=v(2:lx1-1,2:lx2-1,1)+dt*(-w(2:lx1-1,2:lx2-1,1)-v(2:lx1-1,2:lx2-1,1).*(v(2:lx1-1,2:lx2-1,1)-a).*(v(2:lx1-1,2:lx2-1,1)-1)+(D/(dx^2))*(v(3:lx1,2:lx2-1,1)+v(1:lx1-2,2:lx2-1,1)+v(2:lx1-1,3:lx2,1)+v(2:lx1-1,1:lx2-2,1)-4*v(2:lx1-1,2:lx2-1,1)));
     
    % because of periodic BC in the x1 coord, the space averaging must be explicitly defined at near the endpoints for v
    v(1,2:lx2-1,2)=v(1,2:lx2-1,1)+dt*(-w(1,2:lx2-1,1)-v(1,2:lx2-1,1).*(v(1,2:lx2-1,1)-a).*(v(1,2:lx2-1,1)-1)+(D/(dx^2))*(v(2,2:lx2-1,1)+v(lx1-1,2:lx2-1,1)+v(1,3:lx2,1)+v(1,1:lx2-2,1)-4*v(1,2:lx2-1,1))+I(1,2:lx2-1));
    v(lx1,:,2)=v(1,:,2);
    v(lx1-1,2:lx2-1,2)=v(lx1-1,2:lx2-1,1)+dt*(-w(lx1-1,2:lx2-1,1)-v(lx1-1,2:lx2-1,1).*(v(lx1-1,2:lx2-1,1)-a).*(v(lx1-1,2:lx2-1,1)-1)+(D/(dx^2))*(v(lx1,2:lx2-1,1)+v(lx1-2,2:lx2-1,1)+v(lx1-1,3:lx2,1)+v(lx1-1,1:lx2-2,1)-4*v(lx1-1,2:lx2-1,1)));
    % 1st order no-flux boundary conditions at ends of x2 
    v(:,1,2)=v(:,2,2);
    v(:,lx2,2)=v(:,lx2-1,2);
    
    w(1:lx1-1,:,2)=w(1:lx1-1,:,1)+dt*e*(b*v(1:lx1-1,:,1)-g*w(1:lx1-1,:,1)-d);
    
    % BC (periodic at ends of x2, no-flux at ends of x1) for w
    w(lx1,:,2)=w(1,:,2); 
    w(:,1,2)=w(:,2,2);
    w(:,lx2,2)=w(:,lx2-1,2);
    
   
    for x1=1:lx1-1;
        for x2=1:lx2;
          %this sets all sink positions back to zero.
          if p(x1,x2)==0
            v(x1,x2,2)=0;
            w(x1,x2,2)=0;
          end
        end
    end
    
    if t==999;    %more snapshots to store
        v(:,:,4)=v(:,:,2);
        w(:,:,4)=w(:,:,2);
    elseif t==1999;
        v(:,:,5)=v(:,:,2);
        w(:,:,5)=w(:,:,2);
    end

        
    v(:,:,1)=v(:,:,2); %this shifts updated time matriz to previous time matriz spot for next iteration
    w(:,:,1)=w(:,:,2);
end %end of time loop
    

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
