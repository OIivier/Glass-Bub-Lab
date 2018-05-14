function [vpm,vps,vp,t1,t2,v,w]=FHN2d_sinks(simtime,D,a,e,q)
tic
% solves 2D FHN eqns. with periodic boundary conditions using Euler method.
%FHN looks like:    dv/dt=-w-v*(v-a)*(v-1)+D*(d^2/dx^2)v+I
%                   dw/dt=e*(b*v-g*w-d)
%with p heterogeneous cells (the heterogeneity specified in the code).

dx=0.02; % space step (cm) 
dt=0.1; % time step (ms)
% NB: D,dt,dx must be chosen such that 1-4*(dt*D/(dx^2))> = 0 to ensure
% proper averaging.


time=0:dt:simtime;
lt=length(time);
lx=201; %number of cells along each column and row. 



ltI=100;  %how many time steps you wanna apply current for, here 10ms
I=zeros(lx,lx);  % applied current
    I(:,2)=10; 
    
v=zeros(lx,lx,3); %these also define initial conditions
w=zeros(lx,lx,3);
t1=zeros(lx,1);
t2=zeros(lx,1);

p=randmat_best(lx,lx,q);

%D=0.0007; % diffusion coeff. in units of cm^2/ms
%a=0.02;
%e=0.01;
b=0.5;
g=1;
d=0;

for t=1:lt;
    
    
    %define the end of the current injection
    if t==ltI
       I(:,2)=0;
    end
    
    v(2:lx-1,2:lx-1,2)=v(2:lx-1,2:lx-1,1)+dt*(-w(2:lx-1,2:lx-1,1)-v(2:lx-1,2:lx-1,1).*(v(2:lx-1,2:lx-1,1)-a).*(v(2:lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(3:lx,2:lx-1,1)+v(1:lx-2,2:lx-1,1)+v(2:lx-1,3:lx,1)+v(2:lx-1,1:lx-2,1)-4*v(2:lx-1,2:lx-1,1))+I(2:lx-1,2:lx-1));
     
    % because of periodic BC in the x1 coord, the space averaging must be explicitly defined at near the endpoints for v
    v(1,2:lx-1,2)=v(1,2:lx-1,1)+dt*(-w(1,2:lx-1,1)-v(1,2:lx-1,1).*(v(1,2:lx-1,1)-a).*(v(1,2:lx-1,1)-1)+(D/(dx^2))*(v(2,2:lx-1,1)+v(lx-1,2:lx-1,1)+v(1,3:lx,1)+v(1,1:lx-2,1)-4*v(1,2:lx-1,1))+I(1,2:lx-1));
    v(lx,:,2)=v(1,:,2);
    v(lx-1,2:lx-1,2)=v(lx-1,2:lx-1,1)+dt*(-w(lx-1,2:lx-1,1)-v(lx-1,2:lx-1,1).*(v(lx-1,2:lx-1,1)-a).*(v(lx-1,2:lx-1,1)-1)+(D/(dx^2))*(v(lx,2:lx-1,1)+v(lx-2,2:lx-1,1)+v(lx-1,3:lx,1)+v(lx-1,1:lx-2,1)-4*v(lx-1,2:lx-1,1))+I(lx-1,2:lx-1));
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
    
    %rezeroing of sinks & crossing time calculations done in the spatial
    %zone loop
    for i1=2:lx; %start at 2 since periodic BC means 1 row of cells is redundant
        for j1=2:lx; %start at 2 since don't want hetero at stim site
            %rezeroing done here
            if p(i1,j1)==0
                v(i1,j1,2)=0; 
                w(i1,j1,2)=0;
            end

            %crossing times calculated here
            if v(i1,50,2)>0.1 && t1(i1)==0
                t1(i1)=time(t);
            end
            if v(i1,200,2)>0.1 && t2(i1)==0
                t2(i1)=time(t);
            end 
            
        end 
    end
    
    
    
    % want to save the end of pulse for viewing
    if t==round(lt/2)
    v(:,:,3)=v(:,:,2); 
    end 
    
    
    %shifting variables from new in old iterate, to old in mew iterate.
    v(:,:,1)=v(:,:,2); 
    w(:,:,1)=w(:,:,2);
    
end

%median based vp calculator. not that good. 
%if median(t2)==0;
%    vp=0;  
%else delta_t=median(t2)-median(t1);

%mean-based vp calculator
for i1=2:lx
    if t1(i1)==0 || t2(i1)==0
        vp(i1)=0;
    else
        delta_t(i1)=t1(i1)-t2(i1);
        vp(i1)=150*dx/(delta_t(i1)*1e-3); %in cm/s
    end
end

vpm=mean(vp(2:end));
vps=std(vp(2:end));


 toc