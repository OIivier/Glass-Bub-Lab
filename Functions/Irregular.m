%% Load data

% [Data,Lx,Ly] = daRead12('as_3.2p.da', 0);
% [Data,Lx,Ly] = daRead12('MJ4_39.da', 0);
Data = OpenNarayan;

%% Bandpass

[Bandpassed] = bandpass(Data);


%% Plot voltage over time

DataToPlot = Bandpassed;
Fps=40;
iPlot = 1;
StepSize = 15;
figure
xAxis=[1:size(DataToPlot,3)]*(1/Fps);
% x=30;
% y=30;
% xMin = floor(StepSize/2);
% xMax = size(DataToPlot,1);
% yMin = floor(StepSize/2);
% yMax = size(DataToPlot,2);
xMin = 10;
xMax = 70;
yMin=10;
yMax=70;
for x=xMin:StepSize:xMax
    for y=yMin:StepSize:yMax
        subplot(floor(size(DataToPlot,1)/StepSize),floor(size(DataToPlot,2)/StepSize),iPlot)
%         plot(xAxis,AutocorHM( squeeze(DataToPlot(x,y,:)) ) )
        plot(xAxis,squeeze(DataToPlot(x,y,:)) )
        title(['Px(' num2str(x) ',' num2str(y) ')'])
        xlabel('Time (s)')
        xlim([1 120])
        iPlot = iPlot+1;
    end
end

%% Plot FT

DataToPlot = Bandpassed;
iPlot = 1;
StepSize = 20;
figure
Fps=40; %40 Hz
xMin = floor(StepSize/2);
xMax = size(DataToPlot,1);
yMin = floor(StepSize/2);
yMax = size(DataToPlot,2);
for x=xMin:StepSize:xMax
    for y=yMin:StepSize:yMax
        %Do Fourier Transform
        DataFT = abs(fft(squeeze(DataToPlot(x,y,:))));
%         DataFT = abs(fft(AutocorHM(squeeze(DataToPlot(x,y,:)))));
        DataFT = DataFT(1:size(DataToPlot,3)/2); %Discard Half of Points
        FrequencyArray = Fps*(0:size(DataToPlot,3)/2-1)/size(DataToPlot,3);
        subplot(floor(size(DataToPlot,1)/StepSize),floor(size(DataToPlot,2)/StepSize),iPlot)
        plot(FrequencyArray,DataFT )
        title(['Px(' num2str(x) ',' num2str(y) ')'])
%         xlim([0 5])
        xlim([0 Fps/2])
        xlabel('Freq. (Hz)')
        ylabel('Ampl.')
        iPlot = iPlot+1;
    end
end