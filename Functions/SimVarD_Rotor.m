function m = SimVarD_Rotor(nameWithoutExtension, ObstacleX, ObstacleY, ObstacleRadius, ObstacleTime, PacemakerX, PacemakerY, NumberOfBeats);

%% #1 - COMPILE, SET PARAMETERS FOR SIMULATION
cd('/home/vincent/Documents/Glass-Bub-Lab/Simulations/');
mex ../px_ch_VarD_Rotor.cpp 

%%

% nameWithoutExtension = 'Sim_01-Mar-2018_D5_v500_R1_1';
Extension = '.mat';
nameToLoad = strcat(nameWithoutExtension,Extension);
load(nameToLoad, 'Map');
if ~exist('Map','var')
    Map = ReadMap(nameWithoutExtension)
end
% 
% Heatmap3(int16(Map),[0,0,0;1,1,1],0,1);

%%
% 2 - DRAW MAP


NoiseSignal = 0;
NoiseRecord = 0;

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

if exist('PacemakerX','var')
    pdx = PacemakerY-round(size(Map,1)/2); % flipped X and Y for coordinate as we want them
else
%     pdx = ceil(3*size(Map,1)/8);       % horizontal offset from center
    pdx = 10;
end
if exist('PacemakerY', 'var')
    pdy = PacemakerX-round(size(Map,2)/2); % flipped X and Y for coordinate as we want them
else
    pdy = 10;       % vertical offset from center
end


if exist('NumberOfBeats', 'var')
    nBeats = NumberOfBeats;
else
    nBeats = 10;    % number of beats
end


for iMap = 1:(size(Map,1)*size(Map,2))
    iMapX = floor((iMap-1)/size(Map,2))+1;
    iMapY = mod(iMap-1,size(Map,1))+1;
end




% PacemakerCenter=[floor(size(Map,1)/2)+pdx,floor(size(Map,2)/2)+pdy];
% for xx=PacemakerCenter(1)-4:PacemakerCenter(1)+4
%     for yy=PacemakerCenter(2)-4:PacemakerCenter(2)+4
%         if floor(sqrt((xx-PacemakerCenter(1))^2+(yy-PacemakerCenter(2))^2))<=4
%             Map(xx,yy) = 0.5;
%         end
%     end
% end


%%

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

%%
Heatmap3(int16(Map'),gray,0,1);
ind = find(Map~=0);
[PosX,PosY] = ind2sub(size(Map), ind);

%%
VarD = zeros(size(PosX,1),1);
for i=1:size(PosX,1)
    VarD(i)=Map(PosX(i),PosY(i));
end

%%
MapMod=Map;
%

%
if exist('ObstacleX', 'var') && exist('ObstacleY', 'var') && exist('ObstacleRadius', 'var') && ObstacleX~=0 && ObstacleY~=0
    MapMod = addObstacle(MapMod, ObstacleX, ObstacleY, ObstacleRadius);
end

% MapMod = addObstacle(MapMod, 30,20,4,1);
% MapMod = addObstacle(MapMod, 30,60,5,1);
% MapMod = addObstacle(MapMod, 50,40,6,1);
% MapMod = addLine(MapMod,40,5,1);
indMod = find(MapMod~=0);
[PosXMod,PosYMod] = ind2sub(size(MapMod), indMod);

VarDMod = zeros(size(PosXMod,1),1);
for i=1:size(PosXMod,1)
    VarDMod(i)=MapMod(PosXMod(i),PosYMod(i));
end

if exist('ObstacleTime', 'var') && ObstacleTime~=0
    MapModTime = ObstacleTime;
else
    MapModTime = 50;
end

% 2D case is convergent if 2*D*dt/(dx^2) < 1/2;
if 2*D*dt/dx/dx >= 0.5
    error(['Diverging: D*dt/dx/dx = ' num2str(2*D*dt/dx/dx) ' >= 0.5'])
end
%% #4 - RUN MEX CODE
% ======================================================================= %


tic
[m, Toff] = px_ch_VarD_Rotor(Nt, dx, dt, D, wH, wLP, pon, poff, pdx, pdy, nBeats, size(PosX,1), PosX, PosY, NoiseSignal, NoiseRecord, size(PosXMod,1), PosXMod, PosYMod, MapModTime, VarD, VarDMod);
toc


% ======================================================================= %

%Name file
file.name = ['S_'];

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
    toAppend = ['NoB',num2str(NumberOfBeats),'_'];
    file.name = strcat(file.name,toAppend);
end

file.name = strcat(file.name,nameWithoutExtension);

% Save variable
% save(nameWithoutExtension,'Map','MapMod','-v7.3');


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


title_str = [' '];
mechanisms = 'TransientRotor'; 

ftime = 29; % [s]. Final time of the movie.

data = Sensor(m);
data.SampFreq = 1/(dt * tau);
data.MaxActivityAtSpecificPx = max(m(35,69,:));
data.MinActivityAtSpecificPx = min(m(35,69,:));

figureFile = strcat('../FiguresAndGifs/VarD_',file.name);

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

