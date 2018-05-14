%% SpectrHomemade Function
function addSpectro(Signal,Fps,WindowSecond,MarkerFace,MarkerColor)

Length= size(Signal,2); %Number of frame in the signal
WindowFrame=ceil(Fps*WindowSecond); %Size of window, in frames
Overlap=0; %Size of overlap, in seconds (not yet implemented)
nWindows=floor(Length/WindowFrame); %Number of windows
for iWindow=1:nWindows
    Y1 = fft(Signal(1+(iWindow-1)*WindowFrame:iWindow*WindowFrame)); %FFT the signal in that window
    FT = abs(Y1/WindowSecond); %Two sided spectrum of the FT (scale by length of window... not sure necessary)
    FT = FT(1:WindowFrame/2+1); %Single sided spectrum
    FT(2:end-1) = 2*FT(2:end-1); %To "fold back" the second part of the spectrum
    MaxAmplitude = max(FT(2:end));
    for iFreq=1:size(FT,2)
        hold on
        if (iFreq>1) && (FT(iFreq)>MaxAmplitude/2)
            if FT(iFreq)==MaxAmplitude
                plot(iWindow*WindowSecond,((iFreq-1)/size(FT,2))*(Fps/2),MarkerFace,'MarkerSize',exp(2*FT(iFreq)/MaxAmplitude)-1, 'MarkerFaceColor', MarkerColor+(1-sqrt(FT(iFreq)/MaxAmplitude))*([1 1 1]-MarkerColor), 'MarkerEdgeColor',MarkerColor+(1-sqrt(FT(iFreq)/MaxAmplitude))*([1 1 1]-MarkerColor));
            else
                plot(iWindow*WindowSecond,((iFreq-1)/size(FT,2))*(Fps/2),MarkerFace,'MarkerSize',(exp(2*FT(iFreq)/MaxAmplitude)-1)/3, 'MarkerFaceColor', MarkerColor+(1-sqrt(FT(iFreq)/MaxAmplitude))*([1 1 1]-MarkerColor), 'MarkerEdgeColor',MarkerColor+(1-sqrt(FT(iFreq)/MaxAmplitude))*([1 1 1]-MarkerColor));
            end
        end
    end
    xlabel('Time (t)');
    ylabel('f (Hz)');
end

end
