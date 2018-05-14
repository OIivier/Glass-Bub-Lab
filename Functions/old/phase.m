function [tau,d,FF1]=phase(b,pixX,pixY)

%b is the whole data while pxX/Y specify a pizel to make F map for.
tic
%calculate first crossing of autocorrelation function to define period,
%tau
[lx, ly, lt]=size(b);
tau=zeros(lx,ly); %preallocate
F=zeros(lx,ly,lt);

for j=1:lx;
    for i=1:ly;
        c=b(i,j,:);
        d=xcorr(c,'coeff');
        
        for t=lt:2*lt-1;
            if (d(t)>0) && (d(t+1)<0)
            tt=t-lt;
            break
            end
        end
        tau(i,j)=tt; %in timesteps
        s=round(lt/tau(i,j));
        for k=1:s-1;
            F(i,j,k)=c(k*tau(i,j));
        end
    end
end
FF1=nonzeros(F(pixX,pixY,:));
[o,oo,ooo]=size(FF1);
FF2=reshape(F(pixX, pixY, 1:ooo-1),1,ooo-1);
FF3=reshape(F(pixX, pixY, 2:ooo),1,ooo-1);
plot(FF2, FF3,'.')
toc