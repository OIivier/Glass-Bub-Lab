function [c,w]=homoFHN1r(circ,time);

% solves homogeneous 1D FHN eqns.  using Euler method.
%FHN looks like:    dc/dt=-w-c*(c-a)*(c-1)+D*(d^2/dx^2)c+I
%                   dw/dt=e*(b*c-g*w-d)
% time is length of time

%all of the position stuff doesn't matter; i just copied this from FHN1r

D=0; % diffusion coeff.
a=0.139;
e=0.008;
b=1;
g=2.54;
d=0;

dx=2; % space step
dt=1; % time step

t=1:(time/dt)+1;
lt=length(t);
x=1:(circ/dx)+1;
lx=length(x);

I=zeros(lx,lt);  % applied perturbation 
I(2,1)=1; %posn. is arbitrary as long as you use the same one in the end

c=zeros(lx,lt);
%c(lx/2-0.5,1)=1;  % initial condition for c if different from zero

w=zeros(lx,lt);
%w(1:lx,1)=1;   % initial condition for w if different from zero

for t=1:time/dt;
   
    c(2:lx-1,t+1)=c(2:lx-1,t)+dt*(-w(2:lx-1,t)-c(2:lx-1,t).*(c(2:lx-1,t)-a).*(c(2:lx-1,t)-1)+(D/(dx^2))*(c(3:lx,t)-2*c(2:lx-1,t)+c(1:lx-2,t))+I(2:lx-1,t));
    
   
    w(lx,t)=w(1,t); 
    w(1:lx-1,t+1)=w(1:lx-1,t)+dt*e*(b*c(1:lx-1,t)-g*w(1:lx-1,t)-d);
    
end

%subplot(2,1,1); plot(t,c(2,:));
%subplot(2,1,2); plot(t,w(2,:));

c=c(2,:);
w=w(2,:);   %these will be fed into FHN1r as initial conditions.
