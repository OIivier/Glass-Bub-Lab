x0=10;
y0=40;
x1=70;
y1=40;

isRecording=0; %set to 0 if it is a simulation
minFrame = 1;
maxFrame = 10000;
Data0 = squeeze(Data(x0,y0,minFrame:maxFrame));
Data1 = squeeze(Data(x1,y1,minFrame:maxFrame));

[Peaks0, Locations0] = findpeaks(Data0);
[Peaks1, Locations1] = findpeaks(Data1);

Threshold = 0.3;

nPeaks0 = 0;
nPeaks1 = 0;
for iPeak = 1:size(Locations0,1)
    if Data0(Locations0(iPeak)) > Threshold*max(abs(Data0))
        if nPeaks0==0
            firstPeak0 = Locations0(iPeak);
        end
        nPeaks0 = nPeaks0+1;
        lastPeak0 = Locations0(iPeak);
    end
end
if nPeaks0 ~= 0
    Frequency0 = nPeaks0/(lastPeak0-firstPeak0);
end

for iPeak = 1:size(Locations1,1)
    if Data1(Locations1(iPeak)) > Threshold*max(abs(Data1))
        if nPeaks1==0
            firstPeak1 = Locations1(iPeak);
        end
        nPeaks1 = nPeaks1+1;
        lastPeak1 = Locations1(iPeak);
    end
end
if nPeaks1 ~= 0
    Frequency1 = nPeaks1/(lastPeak1-firstPeak1);
end

Period0=1/Frequency0;
Period1=1/Frequency1;

Fps = 40;

Time = [1:size(Data0,1)];
xTitle ='Time (frame)';
yTitle = 'v';

if isRecording==1
    Frequency0=Fps*Frequency0;
    Frequency1=Fps*Frequency1;
    Period0 = 1/Frequency0;
    Period1 = 1/Frequency1;
    Time = Time/Fps;
    xTitle = 'Time (s)';
    yTitle = 'Recorded intensity';
end

mFrequency0 = 1000*Frequency0;
mFrequency1 = 1000*Frequency1;

plot(Time,Data0,'linewidth',2);
hold on;
plot(Time,Data1,'linewidth',2);
xlabel(xTitle);
ylabel(yTitle);
legend0Title = strcat('Pixel (', num2str(x0) ,',', num2str(y0) ,')');
legend1Title = strcat('Pixel (', num2str(x1) ,',', num2str(y1) ,')');
leg1 = legend(legend0Title,legend1Title);
set(gca,'fontsize',20)
set(leg1,'FontSize',20);
