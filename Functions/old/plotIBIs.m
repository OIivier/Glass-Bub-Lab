ct1=load('~/Desktop/fortran/pacemakers/fort.20');
ct2=load('~/Desktop/fortran/pacemakers/fort.21');                                                                                          
dt1=diff(ct1); dt2=diff(ct2);

plot(ct2(2:end), dt2,'b*');
hold
plot(ct1(2:end), dt1,'r*');