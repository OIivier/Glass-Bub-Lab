% % Time correlation over one pixel (autocorrelation)

function Output = AutocorHM(RawData,Frameskip)
% RawData = squeeze(Data(40,45,:));

if ~exist('Frameskip','var')
      Frameskip = 1;
end

DataToAnalyze = zeros(round(size(RawData,1)/Frameskip));

for i=1:floor(size(RawData,1)/Frameskip)
    DataToAnalyze(i) = RawData(i*Frameskip);
end

% DataToAnalyze = RawData;
DataToAnalyze = DataToAnalyze - mean(DataToAnalyze);
TimeRangeToExamine = round(size(DataToAnalyze,1));
CorrelationArray=zeros(1,TimeRangeToExamine);
for tau=1:TimeRangeToExamine
    VProduct =0;
    VSquared =0;
    for t=1:TimeRangeToExamine-tau
        VProduct = VProduct + DataToAnalyze(t)*DataToAnalyze(t+tau);
        VSquared = VSquared + DataToAnalyze(t)*DataToAnalyze(t);
    end
    if VSquared~=0
        CorrelationArray(tau) = VProduct/VSquared;
    end
    %tau
end
Output= CorrelationArray';