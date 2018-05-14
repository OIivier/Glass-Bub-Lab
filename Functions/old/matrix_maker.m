function f1=matrix_maker(f)

lx=80; % size of square matrix  you want

f1=zeros(lx,lx,length(f)/(lx^2));

for i=1:length(f)/(lx^2)
    for j=1:lx
    f1(:,j,i)=f(lx*(j-1)+(lx^2)*(i-1)+1:lx*j+(lx^2)*(i-1));   
    end
end