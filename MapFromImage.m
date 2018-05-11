function Map = MapFromImage(inFile,doSave);
%%

cd('/home/vincent/Documents/Research Thesis/Simulations/');
%%
Map=sum(imread(inFile),3);
Map=Map'/max(max(Map));
Heatmap3(Map',gray,0,1);

%%
if doSave==1
    outFile = strcat('M_',inFile(1:end-4),'.mat');
    save(outFile,'Map','-v7.3');
end