%% SpectrHomemade Function
function SpectrHomemade(Signal_1,Signal_2,Fps_1,Fps_2,WindowSecond_1,WindowSecond_2)
hold off;

Length_1= size(Signal_1,2); %Number of frame in the signal
WindowFrame_1=ceil(Fps_1*WindowSecond_1); %Size of window, in frames
Overlap_1=0; %Size of overlap, in seconds (not yet implemented)
nWindows_1=floor(Length_1/WindowFrame_1); %Number of windows
for iWindow=1:nWindows_1
    Y1 = fft(Signal_1(1+(iWindow-1)*WindowFrame_1:iWindow*WindowFrame_1)); %FFT the signal in that window
    FT_1 = abs(Y1/WindowSecond_1); %Two sided spectrum of the FT (scale by length of window... not sure necessary)
    FT_1 = FT_1(1:WindowFrame_1/2+1); %Single sided spectrum
    FT_1(2:end-1) = 2*FT_1(2:end-1); %To "fold back" the second part of the spectrum
    MaxAmplitude_1 = max(FT_1);
    for iFreq=1:size(FT_1,2)
        hold on
        plot(iWindow*WindowSecond_1,((iFreq-1)/size(FT_1,2))*(Fps_1/2),'ko','MarkerSize',5*sqrt(FT_1(iFreq)/MaxAmplitude_1), 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[1-sqrt(FT_1(iFreq)/MaxAmplitude_1) (1-sqrt(FT_1(iFreq)/MaxAmplitude_1))/4 (1-sqrt(FT_1(iFreq)/MaxAmplitude_1))/4]);
    end
    xlabel('Time (t)');
    ylabel('f (Hz)');
end

Length_2= size(Signal_2,2); %Number of frame in the signal
WindowFrame_2=ceil(Fps_2*WindowSecond_2); %Size of window, in frames
Overlap_2=0; %Size of overlap, in seconds (not yet implemented)
nWindows_2=floor(Length_2/WindowFrame_2); %Number of windows
for iWindow=1:nWindows_2
    Y1 = fft(Signal_2(1+(iWindow-1)*WindowFrame_2:iWindow*WindowFrame_2)); %FFT the signal in that window
    FT_2 = abs(Y1/WindowSecond_2); %Two sided spectrum of the FT (scale by length of window... not sure necessary)
    FT_2 = FT_2(1:WindowFrame_2/2+1); %Single sided spectrum
    FT_2(2:end-1) = 2*FT_2(2:end-1); %To "fold back" the second part of the spectrum
    MaxAmplitude_2 = max(FT_2);
    for iFreq=1:size(FT_2,2)
        hold on
        plot(iWindow*WindowSecond_2,((iFreq-1)/size(FT_2,2))*(Fps_2/2),'ks','MarkerSize',5*sqrt(FT_2(iFreq)/MaxAmplitude_2), 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor',[(1-sqrt(FT_2(iFreq)/MaxAmplitude_2))/4 1-sqrt(FT_2(iFreq)/MaxAmplitude_2) (1-sqrt(FT_2(iFreq)/MaxAmplitude_2))/4]);
    end
    xlabel('Time (t)');
    ylabel('f (Hz)');
end

end
