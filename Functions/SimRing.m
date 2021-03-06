function endDistance = SimRing(wL, showMovie, DistanceCM, nBeats, RealLength, dx, RealDuration, dt, Pacemaker1, Pacemaker2, Ip1, Ip2, wLP1, wLP2, Obstacle, ObstacleWidth, ObstacleTime);



%% COMPILE THE C++ CODE
cd('/home/vincent/Documents/Glass-Bub-Lab/Simulations'); %Change to your local path
mex ../Ring.cpp %Note that the ../ means TwoD.cpp is in the parent directory of Sim2D.m


%% SIMULATION PARAMETERS


if ~exist('RealLength','var')
    RealLength = 2,0; % [cm]
end
if ~exist('dx', 'var')
    dx = 0.01;% [cm]
end

Length = round(RealLength/dx); % [pixels]

if ~exist('RealDuration','var')
    RealDuration = 20; % [seconds]
end
if ~exist('dt','var')
    dt = 0.01; % [seconds]
end
MovieFrameSkip = round(1/dt);
% tau = 0.018; %[seconds]/[movie frame]
nFrames = MovieFrameSkip*RealDuration/dt;
% tau = 0.018; %seconds/[movie frame]
% TotalTime = 10; % seconds
% Treal = floor(TotalTime/tau);% Total movie time in [movie frame]
% dt = 0.025;      % [simulation frame]/[movie frame]
% dt = 40*dx^2;
% dt = 0.025;
Nt = nFrames;
D = 0.0003;      % cm^2/[movie frame]

% dtRealTime = dt*tau;
% D_cm = D/tau;    % [cm^2/s]

% Nt = Treal/dt;     % Total frames

if ~exist('wL', 'var')
    wL = 0.54;     % going 0.3 <= wlp <= 0.5. Affects frequency of pacing.
end
wH = 0.4;%0.6;
pon = 0;
poff = Nt;

if DistanceCM==round(RealLength/2)
    Distance = round(DistanceCM/dx);
elseif DistanceCM<round(RealLength/2)
    Distance = ceil(DistanceCM/dx);
elseif DistanceCM>round(RealLength/2)
    Distance = floor(DistanceCM/dx);
else
    Distance = round(DistanceCM/dx);
end
    
%% FIND MAP FROM FILE

if ~exist('showMovie','var')
    showMovie = 1;
end


Pacemaker1 = round(Length/4);
Pacemaker2 = Pacemaker1+Distance;
while(abs(Pacemaker1-Length)<3||abs(Pacemaker2-Length<3))
    Pacemaker1 = Pacemaker1+1;
    Pacemaker2 = Pacemaker2+1;
end

Map = zeros(Length+4,5);
for i=1:Length;
    Map(i+2,3) = 1;
end


   


%% PACEMAKER SETTINGS
Ip0=20;

if ~exist('Pacemaker1', 'var')
    Pacemaker1 = 10;
end

if ~exist('Pacemaker2', 'var')
    Pacemaker2 = Length+1;
end

if ~exist('wLP1', 'var')
    wLP1 = wL;
end

if ~exist('wLP2', 'var')
    wLP2 = wL;
end

if ~exist('Ip1', 'var')
    Ip1 = Ip0;
end

if ~exist('Ip2', 'var')
    Ip2 = Ip0;
end

if ~exist('nBeats', 'var')
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
if exist('Obstacle', 'var') && exist('ObstacleWidth', 'var') && Obstacle~=0 && ObstacleWidth ~= 0
    for i=1:ObstacleWidth
        MapMod(2+Obstacle-round(ObstacleWidth/2)+i,3)=0;
    end
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
% if 2*D*dt/dx/dx >= 0.5
%     error(['Diverging: D*dt/dx/dx = ' num2str(2*D*dt/dx/dx) ' >= 0.5'])
% end

SizeX = size(Map,1);
SizeY = size(Map,2);

%% RUN MEX CODE
% ======================================================================= %

tic
[m, mW, Toff] = Ring(Nt, dx, dt, D, wH, wL, wLP1, wLP2, pon, poff, Pacemaker1, Pacemaker2, nBeats, size(PosX,1), PosX, PosY, size(PosXMod,1), PosXMod, PosYMod, MapModTime, VCoeff, VCoeffMod, SizeX, SizeY, Ip1, Ip2, Length);
toc


% ======================================================================= %


%% Name file

file.name = ['SimRing_'];

file.name = strcat(file.name,'_L',num2str(Length),'_initD',num2str(Distance));

% if exist('Obstacle','var') && ObstacleWidth~=0
%     toAppend = ['Ob',num2str(Obstacle),'_'];
%     file.name = strcat(file.name,toAppend);
% end
% 
% if exist('ObstacleWidth','var') && ObstacleWidth~=0
%     toAppend = ['ObW',num2str(ObstacleWidth),'_'];
%     file.name = strcat(file.name,toAppend);
% end
% 
% if exist('ObstacleTime', 'var') && ObstacleWidth~=0
%     toAppend = ['ObT',num2str(ObstacleTime),'_'];
%     file.name = strcat(file.name,toAppend);
% end
% 
% if exist('Pacemaker','var')
%     toAppend = ['Pm',num2str(Pacemaker),'_'];
%     file.name = strcat(file.name,toAppend);
% end
% 
% if exist('NumberOfBeats','var')
%     toAppend = ['nBeats',num2str(NumberOfBeats),'_'];
%     file.name = strcat(file.name,toAppend);
% end


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


% title_str = [' '];
% ftime = 29; % [s]. Final time of the movie.
% 
% m=permute(m,[2 1 3]);
% 
% data = Sensor(m);
% data.SampFreq = 1/(dt * tau);
% data.MaxActivityAtSpecificPx = 1.7;
% data.MinActivityAtSpecificPx = -1.95;
% 
% % data.MaxActivityAtSpecificPx = max(m(35,69,:));
% % data.MinActivityAtSpecificPx = min(m(35,69,:));
% 
% figureFile = strcat('../FiguresAndGifs/',file.name);
% 
% figure
% % caxis([min(m(35,69,:)) max(m(35,69,:))])
% caxis([-1.95 1.7])
% Movie(data, 'speed', 60, ...
%             'colormap','parula',...
%             'initialframe', 1 , ...
%             'finalframe', ftime /(dt * tau), ...
%             'save', 'true', ...
%             'frmdelay', 1, ...
%             'title', title_str,...
%             'FileName', figureFile);
%MultFreqSim(m,figureFile,Map);


%%

Data = m(4:Length+3,4,1:MovieFrameSkip:end);
DataW = mW(4:Length+3,4,1:MovieFrameSkip:end);
if showMovie==1
    figure('pos',[500 500 220 200]);

    ax = gca;
    ax.NextPlot = 'replaceChildren';

    cd('/home/vincent/Documents/Glass-Bub-Lab/FiguresAndGifs');

    
    minDiam = ceil(size(Data,1)/pi);
    maxDiam =minDiam+30;
    RingMap = zeros(maxDiam+2,maxDiam+2,size(Data,3));
    centerX = round((maxDiam+2)/2);
    centerY = round((maxDiam+2)/2);

    AngleArray = [2*pi/size(Data,1):2*pi/size(Data,1):2*pi];

    for X=1:size(RingMap,1)
        for Y=1:size(RingMap,2)
            if (centerX-X)^2+(centerY-Y)^2<(maxDiam/2)^2 && (centerX-X)^2+(centerY-Y)^2>(minDiam/2)^2
                Angle = mod(atan2(Y-centerY,X-centerX),2*pi);
                [minValue i] = min(abs(AngleArray-Angle));
                RingMap(X,Y,:)=Data(i,:);
            end
        end
    end

    % for Diameter=minDiam:maxDiam
    %     for i=1:size(Data,1)
    % 
    %         Angle = i*2*pi/size(Data,1);
    %         X = round(centerX + (Diameter/2)*cos(Angle)); 
    %         Y = round(centerY + (Diameter/2)*sin(Angle));
    %         RingMap(X,Y,:)=Data(i,:);
    % 
    %     end
    % end


    nFrames = size(Data,3);
    F(nFrames) = struct('cdata',[],'colormap',[]);


    for j = 1:round(nFrames/1)
        s = surf(RingMap(:,:,j)');
        s.EdgeColor = 'none';
        caxis([-1.95 1.7]);
        view([0 90]);




        drawnow

        F(j) = getframe;

        im = frame2im(F(j)); 
        [imind,cm] = rgb2ind(im,128); 

        %Write to the GIF File 
        if j == 1 
            imwrite(imind,cm,file.name,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,file.name,'gif','WriteMode','append'); 
        end 
    end
end


%%
% 
% figure('pos',[100 500 1000 100]);
% 
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% cd('/home/vincent/Documents/Glass-Bub-Lab/FiguresAndGifs');
% 
% DataSmall = m(:,:,1:40:end);
% DataSmall1 = DataSmall(4:Length+3,4,:);
% DataSmall2 = DataSmall1;
% DataSmall = cat(2,DataSmall1,DataSmall2);
% 
% DataSmallW = mW(:,:,1:40:end);
% DataSmallW1 = DataSmallW(4:Length+3,4,:);
% DataSmallW2 = DataSmallW1;
% DataSmallW = cat(2,DataSmallW1,DataSmallW2);
% 
% nFrames = size(DataSmall,3);
% F(nFrames) = struct('cdata',[],'colormap',[]);
% 
% 
% m = squeeze(DataSmall1);
% mW = squeeze(DataSmallW1);
% 
% 
% m1 = squeeze(m(20,:));
% mW1 = squeeze(mW(20,:));
% 
% m2 = squeeze(m(150,:));
% mW2 = squeeze(mW(150,:));
% 
% v = [min(m1):0.01:max(m1)];
% 
% Beta = 0.7;
% Gamma = 0.5;
% 
% w1 = v-v.^3./3;%dv/dt = 0
% w2 = v+Beta;%./(Gamma.*(wH-wL./(1+exp(-4*v))+wL));%dw/dt =0
% 
% 
% for j = 1:round(nFrames/1)
% %     plot(DataSmall(2:Length-2,4,j)');
%     s = surf(DataSmall(:,:,j)');
%     s.EdgeColor = 'none';
%     caxis([-1.95 1.7]);
%     view([0 90]);
%         
%     
%        
%     
%     drawnow
%     
%     F(j) = getframe;
%       
%     im = frame2im(F(j)); 
%     [imind,cm] = rgb2ind(im,128); 
%     
%     %Write to the GIF File 
%     if j == 1 
%         imwrite(imind,cm,file.name,'gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,file.name,'gif','WriteMode','append'); 
%     end 
% end


%%
% 
% figure('pos',[100 500 1000 1000]);
% 
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% cd('/home/vincent/Documents/Glass-Bub-Lab/FiguresAndGifs');
% 
% DataSmall = m(:,:,1:40:end);
% DataSmall1 = DataSmall(4:Length+3,4,:);
% DataSmall2 = DataSmall1;
% DataSmall = cat(2,DataSmall1,DataSmall2);
% 
% DataSmallW = mW(:,:,1:40:end);
% DataSmallW1 = DataSmallW(4:Length+3,4,:);
% DataSmallW2 = DataSmallW1;
% DataSmallW = cat(2,DataSmallW1,DataSmallW2);
% 
% nFrames = size(DataSmall,3);
% F(nFrames) = struct('cdata',[],'colormap',[]);
% 
% 
% m = squeeze(DataSmall1);
% mW = squeeze(DataSmallW1
% 
% 
% m1 = squeeze(m(20,:));
% mW1 = squeeze(mW(20,:));
% 
% m2 = squeeze(m(380,:));
% mW2 = squeeze(mW(380,:));
% 
% v = [min(m1):0.01:max(m1)];
% 
% Beta = 0.7;
% Gamma = 0.5;
% 
% w1 = v-v.^3./3;%dv/dt = 0
% w2 = v+Beta%./(Gamma.*(wH-wL./(1+exp(-4*v))+wL));%dw/dt =0
% 
% 
% for j = round(5*nFrames/7):round(nFrames/1)
% %     plot(DataSmall(2:Length-2,4,j)');
%     subplot(2,2,[1 2]);
%     s = surf(DataSmall(:,:,j)');
%     s.EdgeColor = 'none';
%     caxis([-1.95 1.7]);
%     view([0 90]);
%     
%     
%     
%     subplot(2,2,3)
%     plot(m1,mW1,'lineWidth',2,'Color', [0.9 0.9 0.9]);
%     plot(v,w1,'lineWidth',2,'Color', 'r');
%     hold on;
%     plot(v,w2,'lineWidth',2,'Color', 'g');
%     plot(m1(j:j+10),mW1(j:j+10),'lineWidth',2,'Color', 'b');
% %     leg1 = legend('$\frac{\partial v}{\partial t}=0$','$\frac{\partial w}{\partial t}=0$','$v$ vs. $w$','Location','northwest');
% %     set(leg1,'Interpreter','latex');
% %     set(leg1,'FontSize',27);
%     
%     
%     subplot(2,2,4)
%     plot(m2,mW2,'lineWidth',2,'Color', [0.9 0.9 0.9]);
%     plot(v,w1,'lineWidth',2,'Color', 'r');
%     hold on;
%     plot(v,w2,'lineWidth',2,'Color', 'g');
%     plot(m2(j:j+10),mW2(j:j+10),'lineWidth',2,'Color', 'b');
% %     leg1 = legend('$\frac{\partial v}{\partial t}=0$','$\frac{\partial w}{\partial t}=0$','$v$ vs. $w$','Location','northwest');
% %     set(leg1,'Interpreter','latex');
% %     set(leg1,'FontSize',27);
%     
%     
%     
%     
%     
%     drawnow
%     
%     F(j) = getframe;
%       
%     im = frame2im(F(j)); 
%     [imind,cm] = rgb2ind(im,128); 
%     
%     %Write to the GIF File 
%     if j == 1 
%         imwrite(imind,cm,file.name,'gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,file.name,'gif','WriteMode','append'); 
%     end 
% end
%%

% 
% DataPoint = squeeze(Data(mod(Pacemaker1+3,size(Data,1)),:));
% DataPointW = squeeze(DataW(mod(Pacemaker1+3,size(DataW,1)),:));
% [peaks locs] = findpeaks(DataPoint);
% threshold = 0.5;
% nPeaks = 0;
% SkipNFirstPeaks = 7;
% 
% 
% for i=1:size(peaks,2)
%     if peaks(1,i)>threshold
%         if nPeaks==SkipNFirstPeaks
%             firstPeak=locs(1,i);
%         end
%         lastPeak=locs(1,i);
%         nPeaks = nPeaks+1;
%     end
% end
% 
% nPeaks = nPeaks-SkipNFirstPeaks;
% period = (lastPeak-firstPeak)/nPeaks;
% speed = Length/period;
% meanW = mean(DataPointW);


%% 
% % 

nPeaks=1;
iFrame=1;
endData = squeeze(Data(:,1,end-1));
while(nPeaks==1 || (abs(endData(end)-max(endData))/max(endData)<0.02 && abs(endData(1)-max(endData))/max(endData)<0.02))
    endData = squeeze(Data(:,1,end-1-(iFrame-1)*20));
    [peaksDist locsDist] = findpeaks(endData);
    threshold=0.5;

    for i=1:size(peaksDist,1)
        if peaksDist(i)>threshold
            if nPeaks==1
                firstPeak=locsDist(i);
            end
            lastPeak=locsDist(i);
            nPeaks = nPeaks+1;
        end
    end
    iFrame = iFrame+1;
end


endDistance = dx*(-abs(lastPeak-firstPeak-Length/2)+Length/2);

% plot(squeeze(Data(1,1,:)));