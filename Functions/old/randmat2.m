function p=randmat2(n,m,q);

for r=1:1000; %first loop tries to make n zeros on average,
               %we will want that average, n to be the number we keep, 
               %so when we find it we will break the loop 
    p=rand(m,q);%generate a mxm matrix with n zeros uniformly distributed 


    for i=1:m;
      for j=1:q;
            if p(i,j)<n/(m*q); %the probability of making a zero
            p(i,j)=0;
            else p(i,j)=1;
            end
      end
    end

if m^2-length(find(p))==n, break, end %the loop stops when i get n zeros 
end





    
    