function MapFromFile = ReadMap(nameWithoutExtension)
Extension = '.mat';
nameToLoad = strcat(nameWithoutExtension,Extension);
DataFromSimulation = load(nameToLoad);
DataFromSimulation = DataFromSimulation.m;
NotAlwaysZero = sum(abs(DataFromSimulation),3);
MapFromFile = ones(size(DataFromSimulation,1),size(DataFromSimulation,2));
for x=1:size(DataFromSimulation,1)
    for y=1:size(DataFromSimulation,2)
        if NotAlwaysZero(x,y) == 0
            MapFromFile(x,y) = 0;
        end
    end
end

MapFromFile = MapFromFile'