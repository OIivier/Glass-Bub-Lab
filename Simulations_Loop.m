%% #1 - COMPILE, SET PARAMETERS FOR SIMULATION
cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/');
mex px_ch.cpp 

%%
% MapFromSimulation = ReadMap('Sim_10-Feb-2018_D7_v900_R1_2');

nameWithoutExtension = 'Sim_01-Mar-2018_D5_v500_R1_1';
Extension = '.mat';
nameToLoad = strcat(nameWithoutExtension,Extension);
DataFromSimulation = load(nameToLoad);
MapFromSimulation = DataFromSimulation.Map;
clear DataFromSimulation
% 
% Heatmap3(int16(MapFromSimulation),[0,0,0;1,1,1],0,1);

%%
% 2 - DRAW MAP
iSim = 1;
iToday = 7;
NEach = 7;
NFibr = 1;
NVert = 1;
NRadius = 1;
for iRadius = 1:NRadius
    for iFibr = 1:NFibr
        for iVert = 1:NVert
            for iEach=iToday:NEach
                %%
%                 Fibr = 0.014; %Chance to have "fibrosis" (block), in [0,1]
%                 Vert = 20; %Vert elongation of fibrosis, in [0,infty[
%                 MaxRadius = 5;
                
                Fibr = 0.005;% + 0.002*iFibr; %Chance to have "fibrosis" (block), in [0,1]
                Vert = 5; %Vert elongation of fibrosis, in [0,infty[
                MaxRadius = 1; %+1*iRadius;
                NoiseSignal = 0;
                NoiseRecord = 0;

                MapDimensions = 82;
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
                pdx = ceil(3*MapDimensions/8);       % horizontal offset from center 
                pdy = 0;       % vertical offset from center
                nBeats = 10;    % number of beats

                
                ScaleMap=1/8;
                BlockArray = [0 0 0; 0 0 0; 0 0 0];
                BlockArray(:,:,2) = [ceil(3*MapDimensions/8) ceil(3*MapDimensions/8) ceil(3*MapDimensions/8);ceil(MapDimensions/2) ceil(MapDimensions/2) ceil(MapDimensions/2);ceil(5*MapDimensions/8) ceil(5*MapDimensions/8) ceil(5*MapDimensions/8)];
                BlockArray(:,:,3) = [ceil(3*MapDimensions/8) ceil(MapDimensions/2) ceil(5*MapDimensions/8);ceil(3*MapDimensions/8) ceil(MapDimensions/2) ceil(5*MapDimensions/8);ceil(3*MapDimensions/8) ceil(MapDimensions/2) ceil(5*MapDimensions/8)];
                Map = false(MapDimensions);
                OriginalObstacles = Map;
                for iMap = 1:MapDimensions^2
                    iMapX = floor((iMap-1)/MapDimensions)+1;
                    iMapY = mod(iMap-1,MapDimensions)+1;
                    if ((iMapX-(MapDimensions/2))^2 + (iMapY-(MapDimensions/2))^2)^(0.5)<((MapDimensions/2)-2)
                        if (1-Fibr)>=rand
                            Map(iMapX,iMapY) = 1;
                        else
                            OriginalObstacles(iMapX,iMapY) = 1;
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


                % for Step=1:ceil(Vert*10)
                %     MapToChange = false(MapDimensions);
                %     for iMapX = 2:MapDimensions-1
                %         for iMapY = 2:MapDimensions-1
                %             if Map(iMapX,iMapY)==0
                %                 if rand<0.5
                %                     if Map(iMapX+1,iMapY)==1
                % %                         MapToChange(iMapX+1,iMapY) = ceil(rand*Vert*10);
                %                             MapToChange(iMapX+1,iMapY) = 1;
                %                     end
                %                 end
                %                 if rand<0.5
                %                     if Map(iMapX-1,iMapY)==1 || MapToChange(iMapX,iMapY)==1
                % %                         MapToChange(iMapX-1,iMapY) = ceil(rand*Vert*10);
                %                         MapToChange(iMapX-1,iMapY) = 1;
                %                     end
                %                 end
                %             end
                %         end
                %     end
                %     
                %     for iMapX=1:MapDimensions
                %         for iMapY=1:MapDimensions
                %             if MapToChange(iMapX,iMapY)==1
                %                 Map(iMapX,iMapY) = 0;
                %             end
                %         end
                %     end
                % end


                for iMap = 1:MapDimensions^2
                    iMapX = floor((iMap-1)/MapDimensions)+1;
                    iMapY = mod(iMap-1,MapDimensions)+1;
                    if OriginalObstacles(iMapX,iMapY)==1
                        Radius = round(rand*MaxRadius);
                %         Radius = MaxRadius;
                        if iMapX-Radius>0 && Radius+iMapX<MapDimensions && iMapY-Radius>0 && Radius+iMapY<MapDimensions
                            for iFibrX=iMapX-Radius:iMapX+Radius
                                for iFibrY=iMapY-Radius:iMapY+Radius
                                    if sqrt((Vert+1)*(iFibrX-iMapX)^2+(iFibrY-iMapY)^2)<=Radius
                                        Map(iFibrX,iFibrY) = 0;
                                    end
                                end
                            end
                        end
                    end
                end

                PacemakerCenter=[floor(size(Map,1)/2)+pdy,floor(size(Map,2)/2)+pdx];
                for xx=PacemakerCenter(1)-4:PacemakerCenter(1)+4
                    for yy=PacemakerCenter(2)-4:PacemakerCenter(2)+4
                        if floor(sqrt((xx-PacemakerCenter(1))^2+(yy-PacemakerCenter(2))^2))<=4
                            Map(xx,yy) = 1;
                        end
                    end
                end


                % Map = addObstacle(Map, 30,40,2);
                % Map = addObstacle(Map, 40,20,2);

                
%                 Heatmap3(int16(Map),hsv,-1,1);
%                 Map = ReadMap('Sim_10-Feb-2018_D7_v900_R1_2');
                Map = MapFromSimulation;
                ind = find(Map==1);
                [PosX,PosY] = ind2sub(size(Map), ind);
%                 MapMod=Map;
                MapMod = addObstacle(Map,68,51,6);
                indMod = find(MapMod==1);
                [PosXMod,PosYMod] = ind2sub(size(MapMod), indMod);
                MapModTime = 35000;
                
                % obstacle
                % figure
                % imagesc(Map)
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

                
                % 2D case is convergent if 2*D*dt/(dx^2) < 1/2;
                if 2*D*dt/dx/dx >= 0.5
                    error(['Diverging: D*dt/dx/dx = ' num2str(2*D*dt/dx/dx) ' >= 0.5'])
                end
                %% #4 - RUN MEX CODE
                % ======================================================================= %
                
                
                tic
                [m, Toff] = px_ch(Nt, dx, dt, D, wH, wLP, pon, poff, pdx, pdy, nBeats, size(PosX,1), PosX, PosY, NoiseSignal, NoiseRecord, size(PosXMod,1), PosXMod, PosYMod, MapModTime);
                toc
                
                
                % ======================================================================= %

                % Name file
                % file.name = ['Simul_' num2str(BlockArray(1,1,1)) '-' num2str(BlockArray(1,2,1)) '-' num2str(BlockArray(1,3,1)) '_' num2str(BlockArray(2,1,1)) '-' num2str(BlockArray(2,2,1)) '-' num2str(BlockArray(2,3,1)) '_' num2str(BlockArray(3,1,1)) '-' num2str(BlockArray(3,2,1)) '-' num2str(BlockArray(3,3,1))];
                % file.name = ['Sim_',datestr(now, 'dd-mmm-yyyy'),'_'];
                file.name = ['Simulations/Sim_',datestr(now, 'dd-mmm-yyyy'),'_D',num2str(int16(Fibr*1000)),'_v',num2str(int16(Vert*100)),'_R',num2str(MaxRadius),'_',num2str(iEach)];

                % Save variable

                save(file.name,'m','Map','MapMod','-v7.3');


                %% View movie
%                 skip_s = 1;  % [tau]
%                 skip_fr = round(skip_s/dt);
%                 
%                 
%                 current = 0;
%                 hf = figure;
%                 for q = 1 : skip_fr : size(m,3)
%                    imagesc(m(:,:,q));
%                    caxis([min(m(35,69,:)) max(m(35,69,:))])
%                    axis square
%                    text(0.8,0.9,sprintf('%.1f/%d s',(q*dt*tau/1e-3/1e3),round(Nt*dt*tau/1e-3/1e3)),'FontSize', 18, 'units','normalized')
%                    hold on
%                 
%                    hold off
%                    getframe(gcf)
%                 
%                 end

                % Save movie
            %     title_str = [num2str(nBeats) ' beats'];
            %     mechanisms = 'TransientRotor'; 
            % 
            %     ftime = 29; % [s]. Final time of the movie.
            % 
            %     data = Sensor(m);
            %     data.SampFreq = 1/(dt * tau);
            %     data.MaxActivityAtSpecificPx = max(m(35,69,:));
            %     data.MinActivityAtSpecificPx = min(m(35,69,:));
            % 
            % 
            % 
            %     figure
            %     caxis([min(m(35,69,:)) max(m(35,69,:))])
            %     Movie(data, 'speed', 60, ...
            %                 'pos',   [200 400 250 250], ...
            %                 'colormap','parula',...
            %                 'initialframe', 1 , ...
            %                 'finalframe', ftime /(dt * tau), ...
            %                 'save', 'true', ...
            %                 'frmdelay', 1, ...
            %                 'title', title_str,...
            %                 'FileName', file.name);
                MultFreqSim(m,file.name,Map);
                clear m;
                iSim = iSim+1;
                fprintf('After %d simulation(s), there is %f Gb of memory left.\n', iSim, MemoryUsage);
            end
        end
    end
end


%% Load simulation & make movie


nameWithoutExtension = 'Sim_02-Mar-2018_D5_v500_R1_7';
Extension = '.mat';
nameToLoad = strcat(nameWithoutExtension,Extension);
DataFromSimulation = load(nameToLoad);
DataFromSimulation = DataFromSimulation.m;

dt = 0.025;
tau = 0.018;
Treal = 1666;
Nt = Treal/dt;
skip_s = 1;  % [tau]
skip_fr = round(skip_s/dt);
current = 0;
%hf = figure;
%for q = 1 : skip_fr : size(DataFromSimulation,3)
%   imagesc(DataFromSimulation(:,:,q));
%   caxis([min(DataFromSimulation(40,13,:)) max(DataFromSimulation(40,13,:))])
%   axis square
%   text(0.8,0.9,sprintf('%.1f/%d s',(q*dt*tau/1e-3/1e3),round(Nt*dt*tau/1e-3/1e3)),'FontSize', 18, 'units','normalized')
%   hold on

%   hold off
%   getframe(gcf)

%end

title_str = [' '];
ftime = 29; % [s]. Final time of the movie.

Data = Sensor(DataFromSimulation);
Data.SampFreq = 1/(dt * tau);
Data.MaxActivityAtSpecificPx = max(DataFromSimulation(40,13,:));
Data.MinActivityAtSpecificPx = min(DataFromSimulation(40,13,:));

figure
caxis([min(DataFromSimulation(40,13,:)) max(DataFromSimulation(40,13,:))])
Movie(Data, 'speed', 60, ...
            'pos',   [200 400 250 250], ...
            'colormap','parula',...
            'initialframe', 1 , ...
            'finalframe', ftime /(dt * tau), ...
            'save', 'true', ...
            'frmdelay', 1, ...
            'FileName', nameWithoutExtension);


%% Load simulation and identify rotors

%Interesting files
%Sim_09-Feb-2018_D1_v0_R2_1
%Sim_09-Feb-2018_D2_v0_R2_2
%Sim_09-Feb-2018_D3_v0_R2_1
%Sim_10-Feb-2018_D7_v500_R2_1
%Sim_10-Feb-2018_D7_v900_R1_2
%Sim_10-Feb-2018_D14_v100_R2_2
%Sim_10-Feb-2018_D14_v500_R1_2
%Sim_10-Feb-2018_D14_v500_R3_1
%Sim_10-Feb-2018_D15_v1000_R2_1
%Sim_10-Feb-2018_D15_v1000_R3_4---
%Sim_10-Feb-2018_D20_v1000_R2_1
%Sim_10-Feb-2018_D20_v1000_R3_1
%Sim_10-Feb-2018_D21_v100_R1_2
%Sim_10-Feb-2018_D21_v900_R3_1 !!!
%Sim_10-Feb-2018_D25_v1000_R2_3
%Sim_10-Feb-2018_D25_v1000_R3_3
%Sim_11-Feb-2018_D19_v300_R3_2

%Sim_12-Feb-2018_D18_v2000_R4_1
%Sim_12-Feb-2018_D14_v2000_R4_2
%Sim_12-Feb-2018_D16_v2000_R3_2
%Sim_12-Feb-2018_D18_v2000_R4_2


% nameWithoutExtension = 'Sim_12-Feb-2018_D18_v2000_R4_2';
% Extension = '.mat';
% nameToLoad = strcat(nameWithoutExtension,Extension);
% DataFromSimulation = load(nameToLoad);
% DataFromSimulation = DataFromSimulation.m;

% SimplifiedData = SimplifyData(DataFromSimulation,1,size(DataFromSimulation,3),10);
% MultFreqIdent(DataFromSimulation,file.name);
% MultFreqData(DataFromSimulation, nameWithoutExtension);
%%
% figure
% plot(squeeze(DataFromSimulation(51,28,:)));
% hold on
% plot(squeeze(DataFromSimulation(63,54,:)));
% plot(squeeze(DataFromSimulation(40,45,:)));
% plot(squeeze(DataFromSimulation(15,50,:)));
        
%%
% hold off
% plot(squeeze(DataFromSimulation(45,59,:)));
% hold on
% plot(squeeze(DataFromSimulation(70,20,:)));
% xlabel('Frame')
% ylabel('Simulated voltage')
% hold on
% plot(squeeze(DataFromSimulation(58,57,:)));


%% Phasemap
% PhaseMap = phaseMap_1(SimplifiedData);

%% Movie Phase
% SimplifiedPhaseMap = SimplifyData(PhaseMap,1,size(PhaseMap,3),5);
% movie_maker_phaseMap(SimplifiedPhaseMap,1,size(SimplifiedPhaseMap,3),1,0,0,0,0,1,0,0);
