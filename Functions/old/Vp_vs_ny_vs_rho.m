rho=[20 30 60 80];
lx1=[100 200 300 400 500];

for j=1:4
    for i=1:5;
    [Vpm(i,j), Vperr(i,j),lVp(i,j)]=FHN2d_h_Vp400y(rho(j),0.0007,lx1(i),800);
    end 
end