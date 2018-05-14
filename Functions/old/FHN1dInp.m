clear
time=100;

% solves 1D FHN eqns. with periodic boundary conditions using Euler method.
%FHN looks like:    dc/dt=-w-c*(c-a)*(c-1)+D*(d^2/dx^2)c+I
%                   dw/dt=e*(b*c-g*w-d)
% time is length of time


% NB: D,dt,dx must be chosen such that 1-2*(dt*D/(dx^2))>0

D=0.0007; % diffusion coeff.
a=0.02;
e=0.01;
b=0.5;
g=1;
d=0;


dx=0.02; % space step
dt=0.1; % time step
circ=1; 

t=1:(time/dt)+1;
lt=length(t);

x=1:(circ/dx)+1;
lx=length(x);

I=zeros(lx,lt);  % applied perturbation 
I(2,1:100)=1;

c=zeros(lx,lt);
w=zeros(lx,lt);

%[c(1:201,1),w(1:201,1)]=homoFHN1r(100,200); %initial conditions are made to look like an action potential in space(see Glass &Josephson 1995). this sets up the wave to go backwards in space

for t=1:time/dt;
     
    c(2:lx-1,t+1)=c(2:lx-1,t)+dt*(-w(2:lx-1,t)-c(2:lx-1,t).*(c(2:lx-1,t)-a).*(c(2:lx-1,t)-1)+(D/(dx^2))*(c(3:lx,t)-2*c(2:lx-1,t)+c(1:lx-2,t))+I(2:lx-1,t));
     %c(1,t+1)=c(2,t);
    % because of periodic BC, the space averaging must be explicitly defined at near the endpoints
    %c(1,t+1)=c(1,t)+dt*(-w(1,t)-c(1,t).*(c(1,t)-a).*(c(1,t)-1)+(D/(dx^2))*(c(2,t)-2*c(1,t)+c(lx-1,t))+I(1,t));
    %c(lx,t)=c(1,t);
    %c(lx-1,t+1)=c(lx-1,t)+dt*(-w(lx-1,t)-c(lx-1,t).*(c(lx-1,t)-a).*(c(lx-1,t)-1)+(D/(dx^2))*(c(lx,t)-2*c(lx-1,t)+c(lx-2,t))+I(lx-1,t));
    
    %w(lx,t)=w(1,t); 
    w(1:lx,t+1)=w(1:lx,t)+dt*e*(b*c(1:lx,t)-g*w(1:lx,t)-d);
    %w(1,t+1)=w(2,t);
end

newplot
subplot(2,1,1);plot(c(2,:));%surface(c); %check out in spacetime
hold on 
plot (c(20,:),'r')
subplot(2,1,2); plot(w(2,:));
hold on 
plot (w(20,:),'r')
% to look at waveform in space plot(1:201,c(:,100))