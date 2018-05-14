function MapName = LoadAndAddObstacle(nameWithoutExtension, CenterX, CenterY, Radius, Intensity)

Extension = '.mat';
nameToLoad = strcat(nameWithoutExtension,Extension);
load(nameToLoad, 'Map');
%%
Map = addObstacle(Map,CenterX,CenterY,Radius,Intensity);
%%

cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Simulations/');
outFile = strcat(nameWithoutExtension,'_add_',num2str(CenterX),'_',num2str(CenterY),'_',num2str(Radius),'_',num2str(Intensity),'.mat');
save(outFile,'Map','-v7.3');

Heatmap3(Map',gray,0,1);

MapName = outFile;