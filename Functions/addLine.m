function MapChanged = addLine(MapInFunction, X, Width, Intensity)

if ~exist('Intensity', 'var')
    Intensity=1;
end

MapChanged = MapInFunction;
for x=X-Width+round(Width/2):X+round(Width/2)-1
    for y=1:round(6*size(MapInFunction,2)/7)
        if MapInFunction(x,y)~=0
            MapChanged(y,x) = 1-Intensity;
        end
    end
end