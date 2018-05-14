%% Summary code to run the filtering codes step by step (Min Ju, 9.20.17)

% [SF,Lx,Ly] = daRead12('MJ4_39.da', 0);
% [SF,Lx,Ly] = daRead12('as_2.1.da', 0);
% [SF,Lx,Ly] = daRead12('Nov29_plate1_5_s.da', 0);
% [SF,Lx,Ly] = daRead12('vinoliv_101_grid.da', 0);
[SF,Lx,Ly] = daRead12('chan1_40.da', 0);
% [SF,Lx,Ly] = daRead12('as_5.5p.da', 0);
% DataFromSimulation = load('Sim_11-Jan-2018_01.mat');
% load('Sim_12-Feb-2018_D16_v2000_R3_2.mat','m');
% load('Sim_10-Feb-2018_D25_v1000_R2_3.mat','m');


% 
% [SF2] = bandpass(SF);
% movie_maker_mod(SF2,1,300,1,0,0,0,0,1,0,1);
% 

%%
MultFreqData(SF);

%% Subtract mean (if daRead12 was not used, e.g. for some simulations)

% DataToAnalyze = SF;
% MatrixOfTheMeans = mean(DataToAnalyze,3);
% for y = 1:size(DataToAnalyze,2)
%      for x = 1:size(DataToAnalyze,1)
%          DataToAnalyze(x,y,:) = DataToAnalyze(x,y,:)-MatrixOfTheMeans(x,y);
%      end
% end
% 
% SF = DataToAnalyze; 



%% Identify standard deviation

 % This is the threshold value that is used to zero the boundary regions.
temp = SF;

MeanMat = zeros(size(temp,1), size(temp,2));
StdMat  = zeros(size(temp,1), size(temp,2));
VarMat  = zeros(size(temp,1), size(temp,2));

iniFrame = 1; % The iniFrame and finFrame can be changed as needed.
finFrame = 4800;

for ii = 1 : size(temp,1)
    for jj = 1 : size(temp,2)
        MeanMat(ii,jj) = mean(temp(ii,jj,iniFrame:finFrame));
        StdMat(ii,jj)  = std(temp(ii,jj,iniFrame:finFrame));
        VarMat(ii,jj)  = var(temp(ii,jj,iniFrame:finFrame));
    end
end

figure
imagesc(StdMat);

axis square;
colorbar; % Plots the standard deviation

%% Confirm the regions considered as outside of the monolayer
th = 24;

StdMat_th = StdMat;
StdMat_th(StdMat_th < th) = 0;
StdMat_th(StdMat_th >= th) = 1;

figure

imagesc(StdMat_th); colorbar; axis square;
grid minor;
xlabel 'x'

%% Mask the regions considered as outside of the monolayer

for i = iniFrame : finFrame %Dish(k).TotalFrames
   SF(:,:,i) = temp(:,:,i).* StdMat_th; % makes the outer edges 0
end

%% Space avg
[SF] = space_avg9_mod(SF);


%% Create a single pixel function with everything
% SinglePixel = zeros(size(DataFromSimulation,3),1);
% FrameStep = 1
% for t = 1:size(DataFromSimulation,3)
%     for x = 1:size(DataFromSimulation,1)
%         for y = 1:size(DataFromSimulation,2)
%             SinglePixel(t) = SinglePixel(t) + DataFromSimulation(x,y,t);
%         end
%     end
% end

%% Discretize values

DataToAnalyze = SF;

Discretized = zeros(size(DataToAnalyze,1),size(DataToAnalyze,2),size(DataToAnalyze,3));
Threshold = max(DataToAnalyze(40,75,:))/2;
for t = 1:size(DataToAnalyze,3)
    for x = 1:size(DataToAnalyze,1)
        for y = 1:size(DataToAnalyze,2)
            if DataToAnalyze(x,y,t)>Threshold
                Discretized(x,y,t)=1;
            else
                Discretized(x,y,t)=0;
            end
        end
    end
end

%% Apply the frequency specific bandpass filter


[SF2] = bandpass(SF);
% DataForSpectro = Discretized;
% %percent = powercalc(SF1,SF2);
% %disp(percent)
% figure(1)
% subplot(2,1,1)
% plot(squeeze(DataForSpectro(50,40,:)));
% title('original data')
% subplot(2,1,2)
% plot(squeeze(DataForSpectro(50,40,:)));
% title('after bandpass filter')
% signal = squeeze(DataForSpectro(50,40,:));
% figure(2)
% % s = spectrogram(signal);
% fps =2300;
% spectrogram(signal,2*fps,fps/2,8000,fps,'yaxis');
% % spectrogram(signal,1200,1100,128,40,'yaxis')
% ylim([0,0.02])
% title('spectrogram')
% % %imagesc(abs(s));
% %plot(s);
% %title('spectrogram of bandpassed signal')
% %disp(powercalc(SF2,SF1));

%% Convert the signal data into phase data

[temp] = phaseMap_1(SF2);



%% Moviemaker 

movie_maker_mod(SF2,1,1000,1,0,0,0,0,1,0,1);
% movie_maker_phaseMap(temp,800,1000,1,0,0,0,0,1,0,0);


%% Moviemaker2

title_str = [' '];
mechanisms = 'TransientRotor'; 

ftime = 29; % [s]. Final time of the movie.
% m=DataFromSimulation;
data = Sensor(m);
data.SampFreq = 1/(dt * tau);
data.MaxActivityAtSpecificPx = max(m(35,69,:));
data.MinActivityAtSpecificPx = min(m(35,69,:));

figureFile = strcat('../FiguresAndGifs/Movie');

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
%%

time=[28.97/66640:28.97/66640:28.97]';
min=5.01;
max=8.49;
scale=round(66640/28.97)
plot(time(min*scale:max*scale),squeeze(m(38,48,min*scale:max*scale)),'LineWidth',2)
xlabel('Time (s.)')
ylabel('Voltage (mV)')
hold on
plot(time(min*scale:max*scale),squeeze(m(40,30,min*scale:max*scale)),'--','LineWidth',2)
set(gcf, 'Position', [0, 0, 500, 300])

%% Find rotors

[NewSourceMap] = windingAll_count(temp, 800, 1, 1, 0) 
%