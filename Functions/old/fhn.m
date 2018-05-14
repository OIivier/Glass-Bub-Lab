a=0.02;
p_c=[0.002,0.01,0.015];
b=0.5;
g=1;
d=0;
e=0.01;

dt=0.001;
t=1:200/dt;
lt=length(t);

ltI=10; %duration of pulse

v=zeros(1,lt);
w=zeros(1,lt);

for c=1:3;
    e=p_c(c);
    for t=1:ltI;
    v(t+1)=v(t)+dt*(-w(t)-v(t).*(v(t)-a).*(v(t)-1)+10);
    w(t+1)=w(t)+dt*e*(b*v(t)-g*w(t)-d);
    end

    for t=ltI:lt;
    v(t+1)=v(t)+dt*(-w(t)-v(t).*(v(t)-a).*(v(t)-1));
    w(t+1)=w(t)+dt*e*(b*v(t)-g*w(t)-d);
    end
    subplot(3,1,c); plot(v);
    hold on
    plot(w,'r');
    
end
