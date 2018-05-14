function MapName = LoadAndAddLine(nameWithoutExtension, X, Width, Intensity)

Extension = '.mat';
nameToLoad = strcat(nameWithoutExtension,Extension);
load(nameToLoad, 'Map');
%%
Map = addLine(Map, X, Width, Intensity);
%%

cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Simulations/');
outFile = strcat(nameWithoutExtension,'_addL_',num2str(X),'_',num2str(Width),'_',num2str(Intensity),'.mat');
save(outFile,'Map','-v7.3');

Heatmap3(Map',gray,0,1);

MapName = outFile;