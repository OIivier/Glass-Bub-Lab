function Map = MapFromImage(inFile,doSave);
%%

%%
Map=sum(imread(inFile),3);
Map=Map'/max(max(Map));
Heatmap3(Map',gray,0,1);