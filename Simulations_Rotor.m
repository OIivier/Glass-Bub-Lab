%% #1 - COMPILE
cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Simulations/Vincent codes');
mex px_ch_Rotor.cpp 

%% #2 - DRAW MAP

NoiseSignal = 0;
NoiseRecord = 0;

MapDimensions = 82;
ScaleMap=1/8;
BlockArray = [0 0 0; 0 0 0; 0 0 0];
BlockArray(:,:,2) = [ceil(3*MapDimensions/8) ceil(3*MapDimensions/8) ceil(3*MapDimensions/8);ceil(MapDimensions/2) ceil(MapDimensions/2) ceil(MapDimensions/2);ceil(5*MapDimensions/8) ceil(5*MapDimensions/8) ceil(5*MapDimensions/8)];
BlockArray(:,:,3) = [ceil(3*MapDimensions/8) ceil(MapDimensions/2) ceil(5*MapDimensions/8);ceil(3*MapDimensions/8) ceil(MapDimensions/2) ceil(5*MapDimensions/8);ceil(3*MapDimensions/8) ceil(MapDimensions/2) ceil(5*MapDimensions/8)];
Map = false(MapDimensions);
OriginalMap = Map;
Fibr = 0.0; %Chance to have "fibrosis" (block), in [0,1]
Vert = 0.0; %Vert elongation of fibrosis, in [0,1]
MaxRadius = 0;
for iMap = 1:MapDimensions^2
    iMapX = floor((iMap-1)/MapDimensions)+1;
    iMapY = mod(iMap-1,MapDimensions)+1;
    if ((iMapX-(MapDimensions/2))^2 + (iMapY-(MapDimensions/2))^2)^(0.5)<((MapDimensions/2)-2)
        if (1-Fibr)>rand
            Map(iMapX,iMapY) = 1;
        else
            OriginalMap(iMapX,iMapY) = 1;
        end
%         for hBlock = 1:3
%             for vBlock = 1:3
%                 if ceil(((iMapX-BlockArray(hBlock,vBlock,2))^2 + (iMapY-BlockArray(hBlock,vBlock,3))^2)^(0.5))<BlockArray(hBlock,vBlock,1)/2
%                     Map(iMapX,iMapY) = 0;
%                 end
%             end
%         end
    end
end


for Step=1:ceil(Vert*10)
    MapToChange = false(MapDimensions);
    for iMapX = 2:MapDimensions-1
        for iMapY = 2:MapDimensions-1
            if Map(iMapX,iMapY)==0
                if rand<0.5
                    if Map(iMapX+1,iMapY)==1
%                         MapToChange(iMapX+1,iMapY) = ceil(rand*Vert*10);
                            MapToChange(iMapX+1,iMapY) = 1;
                    end
                end
                if rand<0.5
                    if Map(iMapX-1,iMapY)==1 || MapToChange(iMapX,iMapY)==1
%                         MapToChange(iMapX-1,iMapY) = ceil(rand*Vert*10);
                        MapToChange(iMapX-1,iMapY) = 1;
                    end
                end
            end
        end
    end
    
    for iMapX=1:MapDimensions
        for iMapY=1:MapDimensions
            if MapToChange(iMapX,iMapY)==1
                Map(iMapX,iMapY) = 0;
            end
        end
    end
end
 

for iMap = 1:MapDimensions^2
    iMapX = floor((iMap-1)/MapDimensions)+1;
    iMapY = mod(iMap-1,MapDimensions)+1;
    if OriginalMap(iMapX,iMapY)==1
        Radius = round(rand*MaxRadius);
%         Radius = MaxRadius;
        if iMapX-Radius>0 && Radius+iMapX<MapDimensions && iMapY-Radius>0 && Radius+iMapY<MapDimensions
            for iFibrX=iMapX-Radius:iMapX+Radius
                for iFibrY=iMapY-Radius:iMapY+Radius
                    if ceil(sqrt((iFibrX-iMapX)^2+(iFibrY-iMapY)^2))<=Radius
                        Map(iFibrX,iFibrY) = 0;
                    end
                end
            end
        end
    end
end






Heatmap2(double(Map),1) %Show the map

% Map = addObstacle(Map, 30,30,4);
% Map = addObstacle(Map, 50,30,2);
% Map = addObstacle(Map, 30,60,6);

px = Map;
ind = find(px==1);
[ipx,jpx] = ind2sub(size(px), ind);

PosX = jpx;
PosY = ipx;

% obstacle
% figure
% imagesc(px)
% colorbar
% axis square
% xlabel('x (px)')
% ylabel('y (px)')

% Generate mesh of positions -------------------------------------------- %
% nx = 1:80;
% ny = 1:80;
% [PosX,PosY] = meshgrid(nx,ny);
% PosX = PosX(:);
% PosY = PosY(:);
%% #3 - SET PARAMETERS FOR SIMULATION

tau = 0.018;     % s

Treal = 1666;    % real time [tau]. Total simulation time in units of tau.
dt = 0.025;      % [tau]
dx = 0.025;      % [cm]
D = 0.0006;      % [cm^2/tau]
D_cm = D/tau;    % [cm^2/s]

Nt = Treal/dt;     % total time points

wLP = 0.48;     % going 0.3 <= wlp <= 0.5. Affects frequency of pacing.
wH = 0.6;
pon = 0;
poff = Nt; 
pdx=0;
% pdx = ceil(3*MapDimensions/8);       % horizontal offset from center 
pdy = 0;       % vertical offset from center
nBeats = 10;    % number of beats

% 2D case is convergent if 2*D*dt/(dx^2) < 1/2;
if 2*D*dt/dx/dx >= 0.5
    error(['Diverging: D*dt/dx/dx = ' num2str(2*D*dt/dx/dx) ' >= 0.5'])
end
%% #4 - RUN MEX CODE
% ======================================================================= %
tic
[m, Toff] = px_ch_Rotor(Nt, dx, dt, D, wH, wLP, pon, poff, pdx, pdy, nBeats, size(PosX,1), PosX, PosY, NoiseSignal, NoiseRecord);
toc
% ======================================================================= %

%% Name file
% file.name = ['Simul_' num2str(BlockArray(1,1,1)) '-' num2str(BlockArray(1,2,1)) '-' num2str(BlockArray(1,3,1)) '_' num2str(BlockArray(2,1,1)) '-' num2str(BlockArray(2,2,1)) '-' num2str(BlockArray(2,3,1)) '_' num2str(BlockArray(3,1,1)) '-' num2str(BlockArray(3,2,1)) '-' num2str(BlockArray(3,3,1))];
% file.name = ['Sim_',datestr(now, 'dd-mmm-yyyy'),'_'];
file.name = ['Sim_',datestr(now, 'dd-mmm-yyyy'),'_D',num2str(int16(Fibr*100)),'_v',num2str(int16(Vert*100)),'_R',num2str(MaxRadius)];

%% Save variable

save(file.name,'m','-v7.3');


%% View movie
% skip_s = 1;  % [tau]
% skip_fr = round(skip_s/dt);
% 
% 
% current = 0;
% hf = figure;
% for q = 1 : skip_fr : size(m,3)
%    imagesc(m(:,:,q));
%    caxis([min(m(35,69,:)) max(m(35,69,:))])
%    axis square
%    text(0.8,0.9,sprintf('%.1f/%d s',(q*dt*tau/1e-3/1e3),round(Nt*dt*tau/1e-3/1e3)),'FontSize', 18, 'units','normalized')
%    hold on
% 
%    hold off
%    getframe(gcf)
% 
% end

%% Save movie
title_str = [num2str(nBeats) ' beats'];
mechanisms = 'TransientRotor'; 

ftime = 29; % [s]. Final time of the movie.

data = Sensor(m);
data.SampFreq = 1/(dt * tau);
data.MaxActivityAtSpecificPx = max(m(35,69,:));
data.MinActivityAtSpecificPx = min(m(35,69,:));



figure
caxis([min(m(35,69,:)) max(m(35,69,:))])
Movie(data, 'speed', 60, ...
            'pos',   [200 400 250 250], ...
            'colormap','parula',...
            'initialframe', 1 , ...
            'finalframe', ftime /(dt * tau), ...
            'save', 'true', ...
            'frmdelay', 1, ...
            'title', title_str,...
            'FileName', file.name);
%% Functions

function MapChanged = addObstacle(MapInFunction, CenterX, CenterY, Radius)
    MapChanged = MapInFunction;
    for x=CenterX-Radius:CenterX+Radius
        for y=CenterY-floor(sqrt(Radius^2-(CenterX-x)^2)):CenterY+floor(sqrt(Radius^2-(CenterX-x)^2))
            MapChanged(x,y) = 0;
        end
    end
end