
NFibrStep = 4
for iSim=1:2*(NFibrStep^2)
    Fibr = 0.01+(0.01*(mod( ((iSim+mod(iSim,2))/2) -1 ,NFibrStep))); %Chance to have "fibrosis" (block), in [0,1]
    Vert = 1+floor( ( ((iSim+mod(iSim,2))/2) -1)/NFibrStep ); %Vert elongation of fibrosis, in [0,infty[
    MaxRadius = 2;
    fprintf('%d. Fibr = %d%%, Vert = %.2f, MaxRadius = %dpx.\n', iSim, int16(Fibr*100), Vert, MaxRadius);
end
