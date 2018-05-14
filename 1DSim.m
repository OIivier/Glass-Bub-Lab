length = 100;
duration = 1000;
v = zeros(length,duration);
w = zeros(length,duration);

for x=1:length;
    v(x,1) = -1.0367;
    w(x,1) = -0.6656;
end

beta = 0.7;
gamma = 0.5;
epsilon = 0.44;
wh = 0.4;
wl = 0.62;
D = 0.0003;

dx = 1;
dt = 1;


for t=1:duration-1;
    for x=3:length-3;
        
        v(x,t+dt) = dt*( (1/epsilon)*(v(x,t)-(v(x,t)^3)/3-w(x,t)) + D*( (v(x+dx,t)-2*v(x,t)+v(x-dx,t)) )/(dx^2)  )+v(x,t);
        w(x,t+dt) = dt*( epsilon*(v(x,t) + beta - gamma*w(x,t))*((wh-wl)/(1+exp(-4*v(x,t))) +wl ) ) + w(x,t);
        
    end
end

v(29,:)