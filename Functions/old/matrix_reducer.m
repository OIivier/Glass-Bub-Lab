function f1=matrix_reducer(f,p)

%takes out every pth matrix form the q from an m*n*q array

lx=80; % size of square matrix  you want
f1=zeros(lx,lx,length(f)/p);

for i=1:length(f)
    if mod(i,10)==0 
     f1(:,:,i/10)=f(:,:,i);
    end
end

