function p=randmatnocorrection(n,m);

%for r=1:1000; %first loop tries to make n zeros on average,
               %we will want that average, n to be the number we keep, 
               %so when we find it we will break the loop 
    %generate a mxm matrix with n zeros uniformly distributed 
p=ones(201);

    for i=1:m;
        for j=1:m
            ii(i)=rand(1)*m;
            jj(j)=rand(1)*m;
            p(round(ii(i)+1),round(jj(j)+1))=0;
      end
    end
%if m^2-length(find(p))==n, break, end %the loop stops when i get n zeros 
end





    
    