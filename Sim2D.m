function m = Sim2D(MapSource, RotorOrPacemaker, PacemakerX, PacemakerY, NumberOfBeats, ObstacleX, ObstacleY, ObstacleRadius, ObstacleTime);

%% COMPILE THE C++ CODE
cd('/home/vincent/Documents/Research Thesis/Simulations'); %Change to your local path
mex ../TwoD.cpp %Note that the ../ means TwoD.cpp is in the parent directory of Sim2D.m

%% FIND MAP FROM FILE

if strcmp(MapSource(end-2:end),'mat')
    load(nameToLoad, 'Map');
    if ~exist('Map','var')
        fprintf('Could not locate the .mat file. Creating a 82x82 blank map instead.\n');
        Map = ones(82,82);
    end

elseif (strcmp(MapSource(end-2:end),'png')||strcmp(MapSource(end-2:end),'jpg')||strcmp(MapSource(end-2:end),'jpeg') )
    if exist(MapSource, 'file') == 2
        Map = MapFromImage(MapSource,0);
    else
        fprintf('Could not locate the image file. Creating a 82x82 blank map instead.\n');
        Map = ones(82,82);
    end

else
    fprintf('Could not locate the image or .mat file. Creating a 82x82 blank map instead.\n');
    Map = ones(82,82);
end
        

%% CHECK WHETHER WE WANT A ROTOR OR A PACEMAKER

if ~exist('RotorOrPacemaker','var')
    fprintf('Pacemaker or rotor not specified, assumed pacemaker.\n');
    RotorOrPaceMaker = 'P';
    isRotor = 0;

elseif strcmp(lower(RotorOrPacemaker),'pacemaker')||strcmp(lower(RotorOrPacemaker),'pace')||strcmp(lower(RotorOrPacemaker),'p')
    RotorOrPaceMaker = 'P';
    isRotor = 0;
    
elseif strcmp(lower(RotorOrPacemaker),'rotor')||strcmp(lower(RotorOrPacemaker),'rot')||strcmp(lower(RotorOrPacemaker),'r')
    RotorOrPaceMaker = 'R';
    isRotor = 1;

else
    fprintf('Could not understand whether a rotor or a pacemaker was demanded, assumed pacemaker.\n');
    RotorOrPaceMaker = 'P';
    isRotor = 0; 
end
    
%% SIMULATION PARAMETERS


tau = 0.018;     % s

Treal = 1666;    % real time [tau]. Total simulation time in units of tau.
dt = 0.025;      % [tau]
dx = 0.025;      % [cm]
D = 0.0003;      % [cm^2/tau]
D_cm = D/tau;    % [cm^2/s]

Nt = Treal/dt;     % total time points

wLP = 0.48*0.5;     % going 0.3 <= wlp <= 0.5. Affects frequency of pacing.
wH = 0.4;%0.6;
pon = 0;
poff = Nt;

%% PACEMAKER SETTINGS

if exist('PacemakerX','var')
    pdx = round(size(Map,1)/2)-PacemakerY+1; %Must switch X and Y to please the convention, +1 adjustment
else
%     pdx = ceil(3*size(Map,1)/8); % vertical offset from center
    pdx = 10;
end
if exist('PacemakerY', 'var')
    pdy = round(size(Map,2)/2)-PacemakerX+1; %Must switch X and Y to please the convention, +1 adjustment
else
    pdy = 10; % horizontal offset from center
end


if exist('NumberOfBeats', 'var')
    nBeats = NumberOfBeats;
else
    nBeats = 10;    % number of beats
end


%% I just keep that here as I might want to re-use it later... it shows how to iterate through the map using only one loop

% for iMap = 1:(size(Map,1)*size(Map,2))
%     iMapX = floor((iMap-1)/size(Map,2))+1;
%     iMapY = mod(iMap-1,size(Map,1))+1;
% end



%% CREATES A BORDER TO THE MAP - Otherwise the C++ code doesn't work

for i=1:size(Map,1)
    Map(i,1)=0;
    Map(i,2)=0;
    Map(i,size(Map,2))=0;
    Map(i,size(Map,2)-1)=0;
end

for i=1:size(Map,2)
    Map(1,i)=0;
    Map(2,i)=0;
    Map(size(Map,1),i)=0;
    Map(size(Map,1)-1,i)=0;
end

%% CREATE ARRAYS OF X AND Y COORDINATES FOR NON-BLOCK LOCATIONS

% Heatmap3(int16(Map'),gray,0,1); %Just shows what the map looks like -- can comment this line
ind = find(Map~=0);
[PosX,PosY] = ind2sub(size(Map), ind);

%% VCoeff is the map of the "gray" areas, where the medium has variable coefficient -- ignore this for simple simulations

VCoeff = zeros(size(PosX,1),1);
for i=1:size(PosX,1)
    VCoeff(i)=Map(PosX(i),PosY(i));
end

%% MapMod is the map after having added an obstacle DURING the simulation. By default, MapMod = Map so the net result is no added osbtacle -- ignore this for simple simulations
MapMod=Map;
%

%
if exist('ObstacleX', 'var') && exist('ObstacleY', 'var') && exist('ObstacleRadius', 'var') && ObstacleX~=0 && ObstacleY~=0
    MapMod = addObstacle(MapMod, ObstacleX, ObstacleY, ObstacleRadius);
end

% MapMod = addObstacle(MapMod, 30,20,4,1); %Example for adding another round obstacle
% MapMod = addLine(MapMod,40,5,1); %Example for adding a horizontal line obstacle

indMod = find(MapMod~=0);
[PosXMod,PosYMod] = ind2sub(size(MapMod), indMod);

VCoeffMod = zeros(size(PosXMod,1),1);
for i=1:size(PosXMod,1)
    VCoeffMod(i)=MapMod(PosXMod(i),PosYMod(i));
end

if exist('ObstacleTime', 'var') && ObstacleTime~=0
    MapModTime = ObstacleTime;
else
    MapModTime = Nt;
end

% 2D case is convergent if 2*D*dt/(dx^2) < 1/2;
if 2*D*dt/dx/dx >= 0.5
    error(['Diverging: D*dt/dx/dx = ' num2str(2*D*dt/dx/dx) ' >= 0.5'])
end

SizeX = size(Map,1);
SizeY = size(Map,2);

%% RUN MEX CODE
% ======================================================================= %


tic
[m, Toff] = TwoD(Nt, dx, dt, D, wH, wLP, pon, poff, pdx, pdy, nBeats, size(PosX,1), PosX, PosY, size(PosXMod,1), PosXMod, PosYMod, MapModTime, VCoeff, VCoeffMod, isRotor, SizeX, SizeY);
toc


% ======================================================================= %


%% Name file

file.name = ['Sim2D_'];

file.name = strcat(file.name,RotorOrPacemaker);

file.name = strcat(file.name,'_');

if exist('ObstacleX','var') && ObstacleRadius~=0
    toAppend = ['ObX',num2str(ObstacleX),'_'];
    file.name = strcat(file.name,toAppend);
end

if exist('ObstacleY','var') && ObstacleRadius~=0
    toAppend = ['ObY',num2str(ObstacleY),'_'];
    file.name = strcat(file.name,toAppend);
end

if exist('ObstacleRadius','var') && ObstacleRadius~=0
    toAppend = ['ObR',num2str(ObstacleRadius),'_'];
    file.name = strcat(file.name,toAppend);
end

if exist('ObstacleTime', 'var') && ObstacleRadius~=0
    toAppend = ['ObT',num2str(ObstacleTime),'_'];
    file.name = strcat(file.name,toAppend);
end

if exist('PacemakerX','var')
    toAppend = ['PmX',num2str(PacemakerX),'_'];
    file.name = strcat(file.name,toAppend);
end

if exist('PacemakerY','var')
    toAppend = ['PmY',num2str(PacemakerY),'_'];
    file.name = strcat(file.name,toAppend);
end

if exist('NumberOfBeats','var')
    toAppend = ['nBeats',num2str(NumberOfBeats),'_'];
    file.name = strcat(file.name,toAppend);
end

file.name = strcat(file.name,MapSource(1:end-4));

% Save map
% save(file.name,'Map','MapMod','-v7.3');


%% View movie (without saving)

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

%% Save movie (also allows you to view it -- must let it run until the end for save to work)


title_str = [' '];
mechanisms = 'TransientRotor'; 

ftime = 29; % [s]. Final time of the movie.

data = Sensor(m);
data.SampFreq = 1/(dt * tau);
data.MaxActivityAtSpecificPx = 1.7;
data.MinActivityAtSpecificPx = -1.95;

% data.MaxActivityAtSpecificPx = max(m(35,69,:));
% data.MinActivityAtSpecificPx = min(m(35,69,:));

figureFile = strcat('../FiguresAndGifs/',file.name);

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
            'FileName', figureFile);
%MultFreqSim(m,figureFile,Map);

