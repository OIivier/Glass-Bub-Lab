function v=randvec(m);

p=rand(m)

for i=1:m
    
if p(i)>0.017
    p(i)=0;
end
end

v=find(p)


    
    