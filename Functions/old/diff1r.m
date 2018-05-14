function c=diff1r(ic, k, dt, dx, time, circ);

%solves 1D diffusion eqn., Ct=k*Cxx with periodic boundary conditions using
%Euler method.

%size of initial condition at c(1,1) is ic
%k is diffusion coeff.
%time & space steps are dt & dx, resp.
%NB: k,dt,dx must be chosen such that 1-2*(dt*k/(dx^2))>0
%time is length of time
%circ is circumference of ring

t=1:(time/dt)+1;
x=1:(circ/dt)+1;

c=zeros(length(x), length(t));
c(1,1)=ic;

for t=1:time/dt;
    c(length(x),t)=c(1,t); %periodic BC
    c(2:length(x)-1,t+1)=c(2:length(x)-1,t)+(dt*k/(dx^2))*(c(3:length(x),t)-2*c(2:length(x)-1,t)+c(1:length(x)-2,t));
    c(1,t+1)=ic-sum(c(2:length(x)-1,t+1));
end
