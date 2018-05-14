v = [-2.0:0.01:2.0];

Beta = 0.7;
Gamma = 0.5;
wL = 0.4;
wH = 0.4;

w1 = v-v.^3./3;%dv/dt = 0
w2 = v+Beta%./(Gamma.*(wH-wL./(1+exp(-4*v))+wL));%dw/dt =0

load('DataVPoint', 'DataVPoint');
load('DataWPoint', 'DataWPoint');
hold off;
plot(v,w1,'lineWidth',2);
hold on;
plot(v,w2,'lineWidth',2);
min = 1;
max = min+1000;
plot(DataVPoint(min:max),DataWPoint(min:max),'lineWidth',2);
leg1 = legend('$\frac{\partial v}{\partial t}=0$','$\frac{\partial w}{\partial t}=0$','$v$ vs. $w$','Location','northwest');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',27);
