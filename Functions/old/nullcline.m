a=0.02;
e=0.01;
b=0.5;
g=1;
d=0;


v_x=linspace(-1,1,100);
w_y=linspace(-1.5,2,100)
w1=-v_x.*(v_x-a).*(v_x-1);
w2=b*v_x-d;


plot(v_x,w1,'r')
hold on
plot(v_x,w2)

[v,w]= meshgrid(-1.5:.2:2);
v_dot=-v.*(v-a).*(v-1)-w;
w_dot=e*(b*v-g*w-d);
quiver(v,w,v_dot, w_dot,2)
