function MapChanged = addObstacle(MapInFunction, CenterX, CenterY, Radius, Intensity)

if ~exist('Intensity', 'var')
    Intensity=1;
end

MapChanged = MapInFunction;
for x=CenterX-Radius:CenterX+Radius
    for y=CenterY-floor(sqrt(Radius^2-(CenterX-x)^2)):CenterY+floor(sqrt(Radius^2-(CenterX-x)^2))
        MapChanged(y,x) = 1-Intensity;
    end
end