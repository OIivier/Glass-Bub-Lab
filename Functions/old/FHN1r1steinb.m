
function [c,w]=FHN1Dsteinb(time);

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
dt=0.0001; % time step
circ=4; %circumference of ring 

% NB this and dt need to be this way for homoFHM1r to work

t=1:(time/dt)+1;
lt=length(t);

x=1:(circ/dx)+1;
lx=length(x);

I=zeros(lx,lt);  % applied perturbation 
I(2,1:100)=1;

c=zeros(lx,lt);
w=zeros(lx,lt);

%[c(1:201,1),w(1:201,1)]=homoFHN1r(100,200); %initial conditions are made to look like an action potential in space(see Glass &Josephson 1995). this sets up the wave to go backwards in space

for t=1:lt;
     
    c(2:lx-1,t+1)=c(2:lx-1,t)+dt*(-w(2:lx-1,t)-c(2:lx-1,t).*(c(2:lx-1,t)-a).*(c(2:lx-1,t)-1)+(D/(dx^2))*(c(3:lx,t)-2*c(2:lx-1,t)+c(1:lx-2,t))+I(2:lx-1,t));
     
   
    w(1:lx-1,t+1)=w(1:lx-1,t)+dt*e*(b*c(1:lx-1,t)-g*w(1:lx-1,t)-d);
    
end

newplot
subplot(2,1,1); plot(x,c(:,20000));
subplot(2,1,2); plot(x,c(:,30000));%check out in spacetime

% to look at waveform in space plot(1:201,c(:,100))