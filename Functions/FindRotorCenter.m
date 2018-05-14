%% Load data

% [Data,Lx,Ly] = daRead12('as_3.2p.da', 0);
[Data,Lx,Ly] = daRead12('MJ4_39.da', 0);

%% Bandpass

[Bandpassed] = bandpass(Data);

%% Make ProbeList

k=1
ProbeList = zeros(1,2)
ProbeDistance = 10;
for i=0:ProbeDistance:40;
    for j=0:ProbeDistance:40;
        ProbeList(k,1) = 20+i
        ProbeList(k,2) = 20+j
        k = k+1
    end
end

%%
DataToAnalyze = Bandpassed;
StepSize=3;
nProbes = size(ProbeList,1);
VelocityArray = zeros(1,1);
Center = [1 1];
StdArray = zeros(1,1);
iStd = 1;
for x=1:StepSize:size(DataToAnalyze,1)
    for y=1:StepSize:size(DataToAnalyze,2)
        Center = [x y];
        iVelocity=1;
        for iProbe=1:size(ProbeList,1)-1
            for jProbe=iProbe+1:size(ProbeList,1)
                ax = ProbeList(iProbe,1)-Center(1);
                ay = ProbeList(iProbe,2)-Center(2);
                bx = ProbeList(jProbe,1)-Center(1);
                by = ProbeList(jProbe,2)-Center(2);
                aDistance = sqrt(ax^2+ay^2);
                bDistance = sqrt(bx^2+by^2);
                if abs((aDistance/bDistance)-1)<1
                    PhaseDifference = PhaseOf(DataToAnalyze(ProbeList(iProbe,1),ProbeList(iProbe,2),:)) - PhaseOf(DataToAnalyze(ProbeList(jProbe,1),ProbeList(jProbe,2),:));
                    Theta = acos( (ax*bx-ay*by)/( aDistance*bDistance ) );
                    if PhaseDifference ~= 0
                        Velocity = Theta/PhaseDifference;
                        VelocityArray([iVelocity]) = Velocity;
                        iVelocity = iVelocity + 1;
                    end
                end
            end
        end
        StdArray([iStd]) = std(VelocityArray);
        iStd = iStd + 1;
    end
end
[MinStd WhereMinStd] = min(StdArray);
xCenter = StepSize*mod(WhereMinStd,(x-1)/StepSize);
yCenter = StepSize*floor(WhereMinStd/( (x-1)/StepSize) );


%% Convert the signal data into phase data

[PhaseMap] = phaseMap_1(Bandpassed);



%% Moviemaker 

movie_maker_mod(Bandpassed,200,400,1,0,0,0,0,1,0,1);
% movie_maker_phaseMap(PhaseMap,800,1000,1,0,0,0,0,1,0,0);