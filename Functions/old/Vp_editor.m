x=20;
for i=1:length(Vp(x,:)) 
if Vp(x,i)<0 || Vp(x,i)>16
Vp(x,i)=0;
end
end 

Vp1=nonzeros(Vp(x,:))

Vpm(x)=mean(Vp1)
Vperr(x)=std(Vp1)
lVp(x)=length(Vp1)