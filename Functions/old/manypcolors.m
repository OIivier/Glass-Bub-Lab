for i=1:100
   p=squeeze(DataMatrix2(900+5*i,:,:));
   pcolor(p)
   shading interp
   f(i)=getframe
   
end