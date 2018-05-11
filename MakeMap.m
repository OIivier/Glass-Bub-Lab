function MapName = MakeMap(Fibr,Vert,MaxRadius,Intensity);

cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Simulations/');


if ~exist('Fibr', 'var')
    Fibr=0;
end

if ~exist('Vert', 'var')
    Vert=0;
end

if ~exist('MaxRadius', 'var')
    MaxRadius=0;
end

if ~exist('Intensity', 'var')
    Intensity=0;
end

MapDimensions = 82;
Map = zeros(MapDimensions);
OriginalObstacles = Map;
for iMap = 1:MapDimensions^2
    iMapX = floor((iMap-1)/MapDimensions)+1;
    iMapY = mod(iMap-1,MapDimensions)+1;
    if((iMapX-(MapDimensions/2))^2 + (iMapY-(MapDimensions/2))^2)^(0.5)<((MapDimensions/2)-2)
        if (1-Fibr/100)>=rand
            Map(iMapX,iMapY) = 1;
        else
            Map(iMapX,iMapY) = 1-Intensity;
            OriginalObstacles(iMapX,iMapY) = 1;
        end
    end
end



for iMap = 1:MapDimensions^2
    iMapX = floor((iMap-1)/MapDimensions)+1;
    iMapY = mod(iMap-1,MapDimensions)+1;
    if OriginalObstacles(iMapX,iMapY)==1
        Radius = round(rand*MaxRadius);
        if iMapX-Radius>0 && Radius+iMapX<MapDimensions && iMapY-Radius>0 && Radius+iMapY<MapDimensions
            for iFibrX=iMapX-Radius:iMapX+Radius
                for iFibrY=iMapY-Radius:iMapY+Radius
                    if sqrt((Vert+1)*(iFibrX-iMapX)^2+(iFibrY-iMapY)^2)<=Radius && ((iFibrX-(MapDimensions/2))^2 + (iFibrY-(MapDimensions/2))^2)^(0.5)<(MapDimensions/2)-2
                        Map(iFibrX,iFibrY) = 1-Intensity;
                    end
                end
            end
        end
    end
end


% Map = addObstacle(Map, 40,20,2);


Heatmap3(Map',gray,0,1);
nameWithoutNumber = ['M_'];

% Name file
if Fibr~=0
    toAppend = ['D',num2str(int16(Fibr*10)),'p10_'];
    nameWithoutNumber = strcat(nameWithoutNumber,toAppend);
    
    if Vert~=0
        toAppend = ['V',num2str(int16(Vert)),'_'];
        nameWithoutNumber = strcat(nameWithoutNumber,toAppend);
    end
    
    if MaxRadius~=0
        toAppend = ['R',num2str(MaxRadius),'_'];
        nameWithoutNumber = strcat(nameWithoutNumber,toAppend);
    end
end


file.name=strcat(nameWithoutNumber,'0')
NumberOfFile = 0
while exist(strcat(file.name,'.mat'), 'file') == 2
    NumberOfFile = NumberOfFile + 1;
    file.name=strcat(nameWithoutNumber,num2str(NumberOfFile));
end
save(file.name,'Map','-v7.3');
MapName = file.name;