%% Load data from file
% [DataFromFile,Lx,Ly] = daRead12('mj_nocells.da', 0);

%% Load data from simulation (not yet fully functional with data but you can try)
% DataFromSimulation = load('Simul_0-0-3_10-10-10_0-0-9.mat');
% DataFromSimulation = DataFromSimulation.m;
Probe_1 = squeeze(DataFromSimulation(50,40,:))';
Probe_2 = squeeze(DataFromSimulation(20,20,:))';
Probe_3 = squeeze(DataFromSimulation(20,60,:))';
Probe_4 = squeeze(DataFromSimulation(60,20,:))';
Probe_5 = squeeze(DataFromSimulation(50,65,:))';
Probe_6 = squeeze(DataFromSimulation(20,40,:))';

%% Discretize values

% Discretized = zeros(size(DataFromSimulation,1),size(DataFromSimulation,2),size(DataFromSimulation,3));
% Threshold = max(DataFromSimulation(40,75,:))/2;
% for t = 1:size(DataFromSimulation,3)
%     for x = 1:size(DataFromSimulation,1)
%         for y = 1:size(DataFromSimulation,2)
%             if DataFromSimulation(x,y,t)>Threshold
%                 Discretized(x,y,t)=1;
%             else
%                 Discretized(x,y,t)=0;
%             end
%         end
%     end
% end

Probe_1_Discretized = squeeze(Discretized(50,40,:))';
Probe_2_Discretized = squeeze(Discretized(20,20,:))';
Probe_3_Discretized = squeeze(Discretized(20,60,:))';
Probe_4_Discretized = squeeze(Discretized(60,20,:))';
Probe_5_Discretized = squeeze(Discretized(50,65,:))';
Probe_6_Discretized = squeeze(Discretized(20,40,:))';

%% Make dummy test functions
x=[0:0.001:5];
y1=sin(2*pi*50*x.*x);
y2=sin(2*pi*50*(x.^0.5).*x);
y3=sin(2*pi*50*sin(x).*x);

%% Call the spectrogram
addSpectro(Probe_1,2300,1.0,'+',[0,0.6,0]);
addSpectro(Probe_2,2300,1.0,'o',[0,0,0.6]);
addSpectro(Probe_3,2300,1.0,'p',[0.6,0,0]);
addSpectro(Probe_4,2300,1.0,'x',[0.6,0,0.6]);
addSpectro(Probe_5,2300,1.0,'s',[0,0.6,0.6]);
addSpectro(Probe_6,2300,1.0,'d',[0.6,0.6,0]);
% addSpectro(Probe_1_Discretized,2300,1.6,'+',[0,0.6,0]);
% addSpectro(Probe_2_Discretized,2300,1.6,'o',[0,0,0.6]);
% addSpectro(Probe_3_Discretized,2300,1.6,'p',[0.6,0,0]);
% addSpectro(Probe_4_Discretized,2300,1.6,'x',[0.6,0,0.6]);
% addSpectro(Probe_5_Discretized,2300,1.6,'s',[0,0.6,0.6]);
% addSpectro(Probe_6_Discretized,2300,1.6,'d',[0.6,0.6,0]);
ylim([0,15]);

