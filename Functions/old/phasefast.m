t=zeros(1,109);
for i=10:length(ct1)-1;
    IBI2(i)=(ct1(i+1)-ct1(i));
    CI2(i)=(ct2(i+1)-ct1(i));
    t(i)=i
end

plot(CI2, IBI2, '*')

a=[t; CI2; IBI2];