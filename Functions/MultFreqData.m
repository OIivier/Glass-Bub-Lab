
%% Finding peak frequency at probes, from bandpass MinJu code
function rotorboolean = MultFreqData(DataName,fileName);
[SF,Lx,Ly] = daRead12(DataName, 0);
% Data = SF(:,:,1:round(size(SF,3)/2));

Data = SF;

Data = bandpass(Data);
DataName = strcat(DataName,'_bp');

%%
SumOverFrames = sum(abs(Data),3);
Map = ones(size(Data,1),size(Data,2));
for x=1:size(Data,1)
    for y=1:size(Data,2)
        if SumOverFrames(x,y)/size(Data,3) < max(max(SumOverFrames))/2
            Map(x,y) = 0;
        end
    end
end

if ~exist('fileName', 'var')
	fileName = strcat('MultFreq_',DataName);
end
MaximumOfAll = max(max(max(abs(Data))));

%%

nx = size(Data,1);
ny = size(Data,2);
nFrames = size(Data,3);
nProbesX = nx;
nProbesY = ny;

if nProbesX==nx
    dx = 1;
else
    dx = floor(nx/(nProbesX+1));
end

if nProbesY==ny
    dy = 1;
else
    dy = floor(ny/(nProbesY+1));
end

Threshold=0.4;

FrequencyArray = zeros(nProbesX,nProbesY);
nPeaksArray = zeros(nProbesX,nProbesY);
peaksDistArray = zeros(nProbesX,nProbesY);
for x = 1:nProbesX
    for y = 1:nProbesY
%         if max(abs(squeeze(Data(x*dx,y*dy,:)))) > 0.1*MaximumOfAll
            PixelSignal = squeeze(Data(x*dx,y*dy,:));
%             PixelSignalFT = abs(fft(PixelSignal));
%             PixelSignalFT = PixelSignalFT(1:floor(size(PixelSignalFT))/2);
            Frameskip=1;
            
            [Peaks,Locations] = findpeaks(PixelSignal, 'MinPeakProminence',max(abs(PixelSignal))*0.05);
            nPeaks = 0;
            for iPeak = 1:size(Locations,1)
                if PixelSignal(Locations(iPeak)) > Threshold*max(abs(PixelSignal))
                    if nPeaks==0
                        firstPeak = Locations(iPeak);
                    end
                    nPeaks = nPeaks+1;
                    lastPeak = Locations(iPeak);
                end
            end
            if nPeaks ~= 0
                nPeaksArray(x,y) = nPeaks;
                peaksDistArray(x,y) = lastPeak-firstPeak;
                FrequencyArray(x,y) = 40*nPeaks/(Frameskip*(lastPeak-firstPeak));
            end
%         end
    end
end

% FrequencyArray

%% Plot single pixels

hold off
plot(squeeze(Data(40,40,:)))
% hold on
% plot(squeeze(Data(35,45,:)))

%% Deciding if rotor or not
fig=figure('Position', [100, 100, 1000, 800]);

subplot(2,2,1);
Heatmap3(FrequencyArray,jet,0.0,5.0);

subplot(2,2,2);
histogram(FrequencyArray);%, 'BinWidth', 3e-5);
xlabel('Frequency (Hz)')
ylabel('Count')
set(gca, 'YScale', 'log')


subplot(2,2,3);
Heatmap3(int16(Map),[0,0,0;1,1,1],0,1);

figureFile = strcat('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/FiguresAndGifs/',fileName);
print(fig,figureFile,'-dpng')
print(fig,figureFile,'-depsc')

% rotorboolean = 1;
% error = 0.15;
% plot(squeeze(m(4*dx,4*dy,:)));
% disp(Norm_peakf);
% averagef = 0;
% counter= 0;
% for x = 1:nProbesX
%     for y = 1:nProbesY
%         if Norm_peakf(x,y) ~= 0
%             averagef = averagef+Norm_peakf(x,y);
%             counter=counter+1;
%         end
%     end
% end
% averagef = averagef/counter;
% problems = 0;
% for xx = 1:nProbesX
%     for yy = 1:nProbesY
%         if Norm_peakf(xx,yy) ~= 0
%             if abs((Norm_peakf(xx,yy)-averagef))/averagef > error
%                 problems = problems+1;
%                 fprintf('Problem at pixel (%d,%d)\n', xx,yy);
%             end
%         end
%     end
% end
% fprintf('problems = %d\n', problems);
% % fprintf('averagef = '+averagef);
% if isequal(Norm_peakf,zeros(nProbesX,nProbesY))||problems>3
%     rotorboolean = 0;
% end
% RotorPresence = {'No Rotor', 'Yes Rotor'};
% fprintf('%s\n', RotorPresence{1+rotorboolean});
% std(std(Norm_peakf))
% plot(squeeze(m(dx*5,dy*6,:)));