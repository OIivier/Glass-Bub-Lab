function p=randmat_best(m,n,q);
%generates matrix of of ones with random entries zeroed
%m,m, specify dimensions of matriz and q specifies how many ranowm poingt u
%want zeroed. 

p=ones(m,n);

for r=1:q
    ii=1+rand(1)*(m-1);
    jj=1+rand(1)*(n-1);
    ii1=round(ii);
    jj1=round(jj);

    while p(ii1,jj1)==0
        ii=1+rand(1)*(m-1);
        jj=1+rand(1)*(n-1);
        ii1=round(ii);
        jj1=round(jj);
    end
    
    p(ii1,jj1)=0;

end







    
    