%% Load data

% [Data,Lx,Ly] = daRead12('as_3.2p.da', 0);
% [Data,Lx,Ly] = daRead12('MJ4_39.da', 0);
DataFromSimulation = load('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Simulations/Vincent codes/Sim_11-Jan-2018_01.mat');
DataFromSimulation = DataFromSimulation.m;

%%
Data = SimplifyData(DataFromSimulation,1,66640,100);



%% Bandpass

[Bandpassed] = bandpass(Data);

%% Convert the signal data into phase data

% [PhaseMap] = phaseMap_1(Data);
[PhaseMapOld] = phaseMap_1_old(Data);

%% Set probes

DataToAnalyze = Data;
% ProbeGroupsD = 20;
% ProbesTopLeft = zeros(floor(size(DataToAnalyze,1)/ProbeGroupsD)*floor(size(DataToAnalyze,2)/ProbeGroupsD),2);
% for i=1:floor(size(DataToAnalyze,1)/ProbeGroupsD)
%     for j=1:floor(size(DataToAnalyze,2)/ProbeGroupsD)
%         ProbesTopLeft((i-1)*floor(size(DataToAnalyze,2)/ProbeGroupsD)+j,:)=[(i-1)*ProbeGroupsD+floor(ProbeGroupsD/2),(j-1)*ProbeGroupsD+floor(ProbeGroupsD/2)];
%     end
% end

ProbesTopLeft = [40,10;40,20;40,30;10,40;20,40;30,40];
ProbesD = 7;
for k=1:size(ProbesTopLeft,1)
    Probes(k,:,:) = [ProbesTopLeft(k,:);ProbesTopLeft(k,:)+[ProbesD 0];ProbesTopLeft(k,:)+[0 ProbesD];ProbesTopLeft(k,:)+[ProbesD ProbesD]];
end

%% Find theta

k=1;


for k=1:size(ProbesTopLeft,1)
    V(1) = ( PhaseOf(DataToAnalyze(Probes(k,1,1),Probes(k,1,2),:))-PhaseOf(DataToAnalyze(Probes(k,2,1),Probes(k,2,2),:)) )/Distance(Probes(k,1,:),Probes(k,2,:));%Velocity according to first index=PhaseDifference/Distance
    V(2) = ( PhaseOf(DataToAnalyze(Probes(k,1,1),Probes(k,1,2),:))-PhaseOf(DataToAnalyze(Probes(k,3,1),Probes(k,3,2),:)) )/Distance(Probes(k,1,:),Probes(k,3,:));
    V(3) = ( PhaseOf(DataToAnalyze(Probes(k,2,1),Probes(k,2,2),:))-PhaseOf(DataToAnalyze(Probes(k,4,1),Probes(k,4,2),:)) )/Distance(Probes(k,2,:),Probes(k,4,:));
    V(4) = ( PhaseOf(DataToAnalyze(Probes(k,3,1),Probes(k,3,2),:))-PhaseOf(DataToAnalyze(Probes(k,4,1),Probes(k,4,2),:)) )/Distance(Probes(k,3,:),Probes(k,4,:));
%     V(1) = AvgPhaseDiff(PhaseMap, Probes(k,1,:), Probes(k,2,:))/Distance(Probes(k,1,:),Probes(k,2,:));%Velocity according to first index=PhaseDifference/Distance
%     V(2) = AvgPhaseDiff(PhaseMap, Probes(k,1,:), Probes(k,3,:))/Distance(Probes(k,1,:),Probes(k,3,:));
%     V(3) = AvgPhaseDiff(PhaseMap, Probes(k,2,:), Probes(k,4,:))/Distance(Probes(k,2,:),Probes(k,4,:));
%     V(4) = AvgPhaseDiff(PhaseMap, Probes(k,3,:), Probes(k,4,:))/Distance(Probes(k,3,:),Probes(k,4,:));
    Theta(k,1) = mod(-atan2d(V(1),V(2)),360);
    Theta(k,2) = mod(-atan2d(V(1),V(3)),360);
    Theta(k,3) = mod(-atan2d(V(4),V(2)),360);
    Theta(k,4) = mod(-atan2d(V(4),V(3)),360);    
end
Theta

%% Find centers

StepSize=1;
StdArray = zeros(1,80*80);
StdArrayMap = zeros(80,80);
for k=1:size(ProbesTopLeft,1)
    iStd = 1;
    for x=1:StepSize:80
        for y=1:StepSize:80
            for i=1:4
                Beta(i) = mod(atan2d(x-Probes(k,i,1),Probes(k,i,2)-y),360);
                DeltaTheta(i) = mod(Theta(i)-Beta(i),360);
                r=Distance(Probes(k,i,:),[x y]);
                m(i) = -(1/r)*log(1-DeltaTheta(i)/90);
            end
            StdArray(iStd) = StdArray(iStd) + std(m);
            StdArrayMap(x,y) = StdArrayMap(x,y) + std(m);
            iStd = iStd + 1;
        end
    end
end

[MinStd WhereMinStd] = min(StdArray);
xCenter = StepSize*floor(WhereMinStd/(x/StepSize) )+1;
yCenter = StepSize*mod(WhereMinStd,y/StepSize);
if yCenter == 0
    yCenter = 80;
end
if WhereMinStd == 4800;
    xCenter = 80;
end

% Heatmap center
Heatmap2(StdArrayMap,1);



%% Moviemaker 

movie_maker_mod(Data,1,size(Data,3),1,0,0,0,0,1,0,1);
% movie_maker_phaseMap(PhaseMapOld,800,1000,1,0,0,0,0,1,0,0);

%%

[NewSourceMap] = windingAll_count(PhaseMap, 200, 10, 1, 0) 