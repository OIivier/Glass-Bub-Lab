a=0.02;
p_c=linspace(0.005,0.04,10); %the parameter range fot the parameter of interest. 
b=0.5;
g=1;
d=0;
e=0.01;

dt=0.001;
time=0:dt:1000;
lt=length(time);



for c=1:10;
   
    e=p_c(c); %define parameter of interest
   
    %initial conditions
    v=0.2; %this is twice the minimum voltage needed to make AP for the nominal parameter set.
    w=0;
        
   for t=1:lt;
       
        v1=v+dt*(-w-v.*(v-a).*(v-1));
        w1=w+dt*e*(b*v-g*w-d);
    
        if v1<0
            R(c)=t*dt;
            break
        end
    
        v=v1;
        w=w1;
    
    end
end

plot(p_c,R, '*--')