
%% Finding peak frequency at probes, from bandpass MinJu code
% function rotorboolean = MultFreqFourier(Data,Map,fileName)
Data = DataFromSimulation; %Comment-swap w/ first line to run non-function
if ~exist('Map','var')
	Map = Data(:,:,1);
end

if ~exist('fileName', 'var')
	fileName = 'No_Name';
end
FractionToAnalyze = 5/10; %Uses all frame from FractionToAnalyze to the end
Data = Data(:,:,round(size(Data,3)*FractionToAnalyze,0):end);



%%


nx = size(Data,1);
ny = size(Data,2);
nProbesX = nx;
nProbesY = ny;
Frameskip=50;

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

FrequencyArray = zeros(nProbesX,nProbesY);
for x = 1:nProbesX
    for y = 1:nProbesY
        PixelSignal = AutocorHM(squeeze(Data(x*dx,y*dy,:)),Frameskip);
%         PixelSignal = squeeze(Data(x*dx,y*dy,:));
        [Peaks,Locations] = findpeaks(PixelSignal);
        nPeaks = size(Locations,1);
        if nPeaks ~= 0
            FrequencyArray(x,y) = 2221*nPeaks/(Frameskip*(Locations(nPeaks)-Locations(1)));
        end;
    end
end

% FrequencyArray


%% Deciding if rotor or not
fig=figure('Position', [100, 100, 500, 800]);

subplot(2,1,1);
Heatmap3(FrequencyArray,spring,0.0,10.0);

subplot(2,1,2);
histogram(FrequencyArray);%, 'BinWidth', 3e-5);
xlabel('Frequency (Hz)')
ylabel('Count')
set(gca, 'YScale', 'log')
% 
% subplot(2,2,3);
% Heatmap3(int16(Map),hsv,-1,1);
% xlim([0.5 5]);
% ylim([0 20]);
% 
% 
% print(fig,fileName,'-dpng')
% print(fig,fileName,'-depsc')

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