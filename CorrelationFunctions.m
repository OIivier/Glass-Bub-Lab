%% Open the data from file

% [DataFromFile,Lx,Ly] = daRead12('mj_nocells.da', 0);
% [DataFromFile,Lx,Ly] = daRead12('as_2.1.da', 0);
% [DataFromFile,Lx,Ly] = daRead12('vinoliv_101_NoLights.da', 0);
[DataFromFile,Lx,Ly] = daRead12('vinoliv_101_grid.da', 0);
% [DataFromFile,Lx,Ly] = daRead12('Vincent14-11-17_LightsOn.da', 0);
% [DataFromFile,Lx,Ly] = daRead12('mj_3_081117.da', 0);
% [DataFromFile,Lx,Ly] = daRead12('MJ4_39.da', 0);

% DataFromSimulation = load('Simul_0-0-3_10-10-10_0-0-9.mat');
% DataFromSimulation = DataFromSimulation.m;

% cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Data Analysis/Data')
% 
% fileID = fopen('Case 1 QRS Subtracted Basket Data.txt','r');
% formatSpec = '%f';
% DataNarayan = fscanf(fileID,formatSpec);
% NFrames = 4000; % number of frames
% Asz = size(DataNarayan,1);
% DataNarayan = reshape(DataNarayan,NFrames,Asz/NFrames);


%% Apply the frequency specific bandpass filter

[SF2] = bandpass(DataFromFile);

%% Subtract mean if necessary

DataToAnalyze = DataFromFile;
MatrixOfTheMeans = mean(DataToAnalyze,3);
for y = 1:size(DataToAnalyze,2)
     for x = 1:size(DataToAnalyze,1)
         DataToAnalyze(x,y,:) = DataToAnalyze(x,y,:)-MatrixOfTheMeans(x,y);
     end
end

%% Identify standard deviation


DataToAnalyze = SF2;
StdMat  = std(DataToAnalyze,1,3);
figure
imagesc(StdMat);
axis square;
colorbar; % Plots the standard deviation
caxis(gca, [0,100]);

%% Confirm the regions considered as outside of the monolayer

Threshold = 30;
Map = StdMat;
Map(Map < Threshold) = 0;
Map(Map >= Threshold) = 1;

figure

imagesc(Map); colorbar; axis square;
grid minor;
xlabel 'x'



%% Spatial Correlation heat map

DataToAnalyze = DataFromFile;
% DataToAnalyze = DataToAnalyze(:,:,1:4000);
x_0=70;
y_0=50;
VProduct = 0;
CorrelationArray = zeros(size(DataToAnalyze,1),size(DataToAnalyze,2));
CorrelationArrayNormalized = CorrelationArray;
xArray = [1:size(DataToAnalyze,1)];
yArray = [1:size(DataToAnalyze,2)];
VSquared = 0;
% for x_0 = 1:size(DataFromFile,1)
%     for y_0 = 1:size(DataFromFile,2)
        for x = 1:size(DataToAnalyze,1)
%             x=x_0;
%             y=y_0;
            for y = 1:size(DataToAnalyze,2)
                for t = 1:size(DataToAnalyze,3)
                    VProduct = VProduct + ( DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x,y,t) );
                    VSquared = VSquared + abs( DataToAnalyze(x,y,t)*DataToAnalyze(x,y,t) );
                end
                CorrelationArray(x,y)= CorrelationArray(x,y) + VProduct;
                CorrelationArrayNormalized(x,y)= CorrelationArray(x,y)/VSquared;
                VProduct = 0;
                VSquared = 0;
            end
        end
%     end
% end
% VSquared = VSquared/(size(DataToAnalyze,1)*size(DataToAnalyze,2)*size(DataToAnalyze,3));
% CorrelationArrayNormalized = CorrelationArray./VSquared;


%Heatmap 1
% CorrelationArray = flipud(CorrelationArray);
% hm = HeatMap(CorrelationArray);
% close all hidden
% ax = hm.plot;
% colorbar('Peer', ax);
% Scale = 2e3;
% caxis(ax, [-Scale,Scale]);
%%
% Heatmap2(CorrelationArrayNormalized,1,0);
% figure;
% imagesc(CorrelationArrayNormalized,[-1,1]); 
% colormap(GBRColorMap(256)); 
Heatmap3(CorrelationArrayNormalized,GBRColorMap(256),-1.0,1.0);

%% Spatial Correlation All Map

DataToAnalyze = DataFromFile;
DataToAnalyze = SimplifyData(DataToAnalyze,1,size(DataToAnalyze,3),100);
VProduct = 0;
CorrelationArray = zeros(size(DataToAnalyze,1),size(DataToAnalyze,2));
CorrelationArrayNormalized = CorrelationArray;
CorrelationMapTotal = CorrelationArray;
xArray = [1:size(DataToAnalyze,1)];
yArray = [1:size(DataToAnalyze,2)];
VSquared = 0;
for x_0 = 1:size(DataFromFile,1)
    for y_0 = 1:size(DataFromFile,2)
        for x = 1:size(DataToAnalyze,1)
            for y = 1:size(DataToAnalyze,2)
                for t = 1:size(DataToAnalyze,3)
                    VProduct = VProduct + ( DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x,y,t) );
                    VSquared = VSquared + abs( DataToAnalyze(x,y,t)*DataToAnalyze(x,y,t) );
                end
                CorrelationArray(x,y)= CorrelationArray(x,y) + VProduct;
                CorrelationArrayNormalized(x,y)= CorrelationArray(x,y)/VSquared;
                VProduct = 0;
                VSquared = 0;
            end
        end
        CorrelationMapTotal = CorrelationMapTotal + CorrelationArrayNormalized;
        fprintf('(%d,%d)\n', x_0, y_0);
    end
end
CorrelationMapTotalNormalized = CorrelationMapTotal/(size(DataToAnalyze,1)*size(DataToAnalyze,2));

%%
Heatmap2(CorrelationMapTotalNormalized,1,1);



%% Scatter plot

% scatter(squeeze(DataToAnalyze(40,40,:)),squeeze(squeeze(DataToAnalyze(39,41,:))))
% scatter(squeeze(DataToAnalyze(40,40,:)),squeeze(squeeze(DataToAnalyze(40,41,:))))
scatter(squeeze(DataToAnalyze(40,40,:)),squeeze(squeeze(DataToAnalyze(39,41,:))))

%% Time correlation over one pixel (autocorrelation)
DataToAnalyze = DataFromFile;
DataToAnalyze = DataToAnalyze - mean(DataToAnalyze);
TimeRangeToExamine = floor(size(DataToAnalyze,3));
x_0=40;
y_0=40;
CorrelationArray=zeros(1,TimeRangeToExamine);
for tau=1:TimeRangeToExamine
    VProduct =0;
    VSquared =0;
    for t=1:size(DataToAnalyze,3)-tau
        VProduct = VProduct + DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x_0,y_0,t+tau);
        VSquared = VSquared + DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x_0,y_0,t);
    end
    CorrelationArray(tau) = VProduct/VSquared;
    tau
end
TimeAxis=[1/(TimeRangeToExamine):1/(TimeRangeToExamine):1]*30;
plot(TimeAxis,CorrelationArray)
% figure
% plot(TimeAxis,squeeze(DataToAnalyze(40,40,:)))

% For pixel (40,40) in MJ4_39.da, whole range, tau step = 1
% plot(abs(fft(squeeze(DataToAnalyze(40,40,1:1000)))))
% plot(abs(fft(CorrelationArray(1:1000))))


%% Time correlation over one pixel (autocorrelation) Narayan version
DataToAnalyze = DataNarayan(:,11)';
DataToAnalyze = DataToAnalyze - mean(DataToAnalyze);
TimeRangeToExamine = 3000;
CorrelationArray=zeros(1,TimeRangeToExamine);
for tau=1:TimeRangeToExamine
    VProduct =0;
    VSquared =0;
    for t=1:size(DataToAnalyze,2)-tau
        VProduct = VProduct + DataToAnalyze(1,t)*DataToAnalyze(1,t+tau);
        VSquared = VSquared + DataToAnalyze(1,t)*DataToAnalyze(1,t);
    end
    CorrelationArray(tau) = VProduct/VSquared;
    tau
end
plot(CorrelationArray)

% For pixel (40,40) in MJ4_39.da, whole range, tau step = 1
% plot(abs(fft(squeeze(DataToAnalyze(40,40,1:1000)))))
% plot(abs(fft(CorrelationArray(1:1000))))

%% Correlation over a small circle or other geometries

DataToAnalyze = DataFromSimulation;
SpaceRangeToExamine = 8;
x_0=70;
y_0=40;
V_0Avg=mean(DataToAnalyze(x_0,y_0,:));
CorrelationArray=zeros(1,size(DataToAnalyze,3));
for t=1:size(DataToAnalyze,3)
    VProduct=0;
    VSquared=0;
    for x=x_0-SpaceRangeToExamine:x_0+SpaceRangeToExamine
%     for x=x_0:x_0+SpaceRangeToExamine
        for y=y_0-ceil(sqrt(SpaceRangeToExamine^2-(x_0-x)^2)):y_0+ceil(sqrt(SpaceRangeToExamine^2-(x_0-x)^2))
%         for y=y_0:y_0+ceil(sqrt(SpaceRangeToExamine^2-(x_0-x)^2))
%           for y=y_0-1:y_0+1
            if Map(x,y)~=0
                VProduct = VProduct + V_0Avg*DataToAnalyze(x,y,t);
                VSquared = VSquared + DataToAnalyze(x,y,t)*DataToAnalyze(x,y,t);
            end
        end
    end
    CorrelationArray(t) = VProduct/VSquared;
    t
end

TimeAxis=[1/(size(CorrelationArray,2)):1/(size(CorrelationArray,2)):1]*30;
plot(TimeAxis,CorrelationArray)
% figure
% plot(squeeze(DataToAnalyze(40,40,:)))

%% Correlation over a small ring

DataToAnalyze = SF2;
RingInnerRadius = 3;
RingWidth = 2;
x_0=60;
y_0=40;
V_0Avg=mean(DataToAnalyze(x_0,y_0,:));
CorrelationArray=zeros(1,size(DataToAnalyze,3));

for t=1:size(DataToAnalyze,3)
    VProduct=0;
    VSquared=0;
    for RingRadius=RingInnerRadius:RingInnerRadius+RingWidth;
        for x=x_0-RingInnerRadius:x_0+RingRadius
            for y=y_0-ceil(sqrt(RingRadius^2-(x_0-x)^2)):2*ceil(sqrt(RingRadius^2-(x_0-x)^2)):y_0+ceil(sqrt(RingRadius^2-(x_0-x)^2))
                if Map(x,y)~=0
                    VProduct = VProduct + DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x,y,t);
                    VSquared = VSquared + DataToAnalyze(x,y,t)*DataToAnalyze(x,y,t);
                end
            end
        end
    end
    CorrelationArray(t) = VProduct/VSquared;
    t
end

TimeAxis=[1/(size(CorrelationArray,2)):1/(size(CorrelationArray,2)):1]*30;
plot(TimeAxis,CorrelationArray)
% figure
% plot(squeeze(DataToAnalyze(40,40,:)))


%% Correlation over radius of a ring at a fixed pixel in a fixed frame

DataToAnalyze = SF2;
%show frame (65098 is good for pixel (40,28))
%imagesc(squeeze(DataToAnalyze(:,:,300)));
t = 420; %Or any frame that has enough waves
x_0=18;
y_0=40;
Radius = 1;
xSize = size(DataToAnalyze,1);
ySize = size(DataToAnalyze,2);
MaxRadius = min([x_0,y_0,xSize-x_0,ySize-y_0])-1;
V_0Avg=mean(DataToAnalyze(x_0,y_0,:));
CorrelationArray=zeros(1,MaxRadius);

for Radius=1:MaxRadius
    VProduct=0;
    VSquared=0;
    for x=x_0-Radius:x_0+Radius
        for y=y_0-ceil(sqrt(Radius^2-(x_0-x)^2)):y_0+ceil(sqrt(Radius^2-(x_0-x)^2))
            if Map(x,y)~=0
                VProduct = VProduct + DataToAnalyze(x_0,y_0,t)*DataToAnalyze(x,y,t);
                VSquared = VSquared + DataToAnalyze(x,y,t)*DataToAnalyze(x,y,t);
            end
        end
    end
    CorrelationArray(Radius) = VProduct/VSquared;
    Radius
end

RadiusAxis=[1:Radius];
plot(RadiusAxis,CorrelationArray)
% figure
% plot(squeeze(DataToAnalyze(40,40,:)))


%% V Magnitude over a small area

DataToAnalyze = DataFromFile;
SpaceRangeToExamine = 8;
x_0=60;
y_0=40;
V_0Avg=mean(DataToAnalyze(x_0,y_0,:));
MagnitudeOverTimeArray=zeros(1,size(DataToAnalyze,3));
for t=1:size(DataToAnalyze,3)
    VMagnitudeInArea=0;
    SizeOfArea=0;
    for x=x_0-SpaceRangeToExamine:x_0+SpaceRangeToExamine
%     for x=x_0:x_0+SpaceRangeToExamine
        for y=y_0-ceil(sqrt(SpaceRangeToExamine^2-(x_0-x)^2)):y_0+ceil(sqrt(SpaceRangeToExamine^2-(x_0-x)^2))
%         for y=y_0:y_0+ceil(sqrt(SpaceRangeToExamine^2-(x_0-x)^2))
%         for y=y_0-1:y_0+1
            if Map(x,y)~=0
                VMagnitudeInArea = VMagnitudeInArea + DataToAnalyze(x,y,t);
                SizeOfArea = SizeOfArea+1;
            end
        end
    end
    MagnitudeOverTimeArray(t) = VMagnitudeInArea/SizeOfArea;
    t
end

TimeAxis=[1/(size(MagnitudeOverTimeArray,2)):1/(size(MagnitudeOverTimeArray,2)):1]*30;
plot(TimeAxis,MagnitudeOverTimeArray)
% figure
% plot(squeeze(DataToAnalyze(40,40,:)))




%% Moviemaker 

movie_maker_mod(DataFromFile,1000,size(DataFromFile,3),1,0,0,0,0,1,0,1);