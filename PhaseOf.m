function Phase = PhaseOf(DataToAnalyze)
% Find Peak Frequency
DataToAnalyze=squeeze(DataToAnalyze);

%% 1st method

% FT = abs(fft(DataToAnalyze));
% FT = FT(1:size(DataToAnalyze,1)/2);
% [PeakFreqMag, PeakFreq] = max(FT);
% PeakFreq = PeakFreq/4800;
% Period = 1/PeakFreq;
% [NextMax, NextMaxFrame] = max(DataToAnalyze(1:round(Period)));
% Phase = 2*pi*(1-NextMaxFrame/Period);

%% I don't understand this one
% Phase = atan2d(DataToAnalyze(1),DataToAnalyze(1+floor(Period/4)));

%% 2nd method

[Peaks, Locations] = findpeaks(DataToAnalyze);
if size(Peaks,1)==0
    Phase =0;
else
    Phase=-Locations(1);
end