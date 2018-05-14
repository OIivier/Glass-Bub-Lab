function [tau,APD]=risetime(v)

for i=1:length(v(:,1))-1
    if v(i,1)>0.1
       t1=v(i,2);
    break
    end
end

for i=1:length(v(:,1))-1
    if v(i+1,1)<v(i,1) 
        t2=v(i,2);
        break
    end
end


for i=1:length(v(:,1))-1
    if v(i,1)<0
        t3=v(i,2);
        break
    end
end

    
tau=t2-t1;
APD=t3-t1;