function SpatialCorrelation(NameOfRecording,x_0,y_0,Scale)

[Data,Lx,Ly] = daRead12(NameOfRecording, 0);

DataToAnalyze = Data;
% DataToAnalyze = DataToAnalyze(:,:,1:4000);

if ~exist('x_0', 'var')
    x_0=70;
end

if ~exist('y_0', 'var')
    y_0=50;
end

if ~exist('Scale', 'var')
    Scale = 1.0;
end

VProduct = 0;
CorrelationArray = zeros(size(DataToAnalyze,1),size(DataToAnalyze,2));
CorrelationArrayNormalized = CorrelationArray;
xArray = [1:size(DataToAnalyze,1)];
yArray = [1:size(DataToAnalyze,2)];
VSquared_0 = 0;
VSquared = 0;
% for x_0 = 1:size(DataToAnalyze,1)
%     for y_0 = 1:size(DataToAnalyze,2)
for x = 1:size(DataToAnalyze,1)
%     x=x_0;
%     y=y_0;
    for y = 1:size(DataToAnalyze,2)
        for t = 1:size(DataToAnalyze,3)
            VProduct = VProduct + DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x,y,t);
            VSquared_0 = VSquared_0 + DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x_0,y_0,t);
            VSquared = VSquared + DataToAnalyze(x,y,t)*DataToAnalyze(x,y,t);
        end
        CorrelationArray(x,y)= CorrelationArray(x,y) + VProduct;
        CorrelationArrayNormalized(x,y)= CorrelationArray(x,y)/sqrt(VSquared*VSquared_0);
        VProduct = 0;
        VSquared = 0;
        VSquared_0 = 0;
    end
end

Heatmap3(CorrelationArrayNormalized,GBRColorMap(256),-Scale,Scale);