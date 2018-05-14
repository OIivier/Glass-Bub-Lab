function f1=matrix_maker1D(f,lx);

% lx size spatial scale of square matrix you want

f1=zeros(lx,length(f)/(lx));

for i=1:length(f)/lx;
    f1(:,i)=f(lx*(i-1)+1:lx*i);   
end

%f1=(fliplr(f1))';
f1=f1';
pcolor(f1);
shading interp;