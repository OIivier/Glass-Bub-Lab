a=0.02;
e=0.04;
b=0.5;
g=1;
d=0;

dt=0.001;
t=1:500/dt;
lt=length(t);

ltI=10; %duration of pulse

v=zeros(1,lt);
w=zeros(1,lt);


for t=1:ltI;
    v(t+1)=v(t)+dt*(-w(t)-v(t).*(v(t)-a).*(v(t)-1)+20);
    w(t+1)=w(t)+dt*e*(b*v(t)-g*w(t)-d);
end

for t=ltI:lt;
    v(t+1)=v(t)+dt*(-w(t)-v(t).*(v(t)-a).*(v(t)-1));
    w(t+1)=w(t)+dt*e*(b*v(t)-g*w(t)-d);
end
plot(v);
hold on
plot(w,'r');
%plot(v,w)