% Class used to extract information from data file. 
%
% Note: The SF1 matrix is a nRow-by-nCol-by-nPx matrix, which means that SF1(60,
% 20, 100) is the activity of the 60th row, 20th column at the 100th frame.
% This is not the same as pixel (60,20). In fact, this is the activity of
% the pixel (20,60). This version takes this issue into account.
%
%
% L. Campanari, Feb 2016

classdef Sensor
    
    properties (Access = public)
        SampFreq = 40;              % Sampling frequency (Hz)
        StimuliThreshold = 2500;    % Threshold to consider stimulus (a.u.)
    end
    
    properties (Access = public)
        SpecificPxForActivity = [20,60];    % Px used to get the max/min 
                                            % activity, which is used to
                                            % define caxis in 'Movie'. 
                                            % 41,20
                                            % [65,45]
    end
    
    properties (Dependent)
        SampPeriod                  % Sampling period (s)
        TotalTime                   % Total time   recorded (s)
        TotalFrames                 % Total frames recorded
        TimeAxis                    % Time array (s) of full interval
        FLAG_STIM                   % 0: no stimuli / 1: there is stimuli
    end
    
    properties
        % File-related
        Name                        % Dish name
        Date                        % Experiment date
        Treatment                   % Treatment received
        
        % Image-related
        DataMatrix                  % SF1 matrix with corrected axis (SF1 transposed)
        TimeSeries
        FullTimeSeries
        VarMap                      % Variance map of movie (a.u.)
        Monolayer                   % structure
      
        % Stimulation-related
        Stimuli                     % Stimuli time series
        StimuliPeriod
        StimuliFrequency
        StimuliRange
        StimuliReference            % Array with times when depo happens
        StimuliLocation
        
        % Movie-related
        % Used to define the max and min (color) axis for Movie
        MaxActivityAtSpecificPx       % Maximum activity in a specific pixel
        MinActivityAtSpecificPx       % Minimum activity in a specific pixel
        
        % Analysis-related
        xCorr
        lag
        APThreshold                 % Threshold for action potential.
        APCrossTimes                % Time when crossing APThreshold.
        IBI                         % Cell array (eg IBI{3}(:,1) and 
                                    % IBI{3}(:,2) are the IBI values and
                                    % IBI times of coord 3, respectively
    end
    
    methods
        %% functions that get proprieties
        function ThisObj = Sensor(SF1, Stim)
            
%             FixedSF1 = permute(SF1,[2 1 3]); % Transpose SF1
              FixedSF1 = SF1;
            
            if nargin == 1
                ThisObj.DataMatrix = FixedSF1;
            elseif nargin == 2
                ThisObj.DataMatrix = FixedSF1;
                
                ThisObj.Stimuli = Stim;
                [ThisObj.StimuliPeriod, ...
                ThisObj.StimuliRange]   = GetStimInfo(ThisObj);
                ThisObj.StimuliFrequency = 1./ThisObj.StimuliPeriod;
                ThisObj.StimuliReference = ComputeStimReference(ThisObj);
                
                if ThisObj.FLAG_STIM == 0
                    ThisObj.Stimuli = 0;
                    ThisObj.Stimuli(1) = [];
                end
            end
            
            % Get values for max and min to use in caxis (on the movie).
            xk = ThisObj.SpecificPxForActivity(1);
            yk = ThisObj.SpecificPxForActivity(2);
            try
                ThisObj.MaxActivityAtSpecificPx = max(ThisObj.DataMatrix(yk,xk,:));
                ThisObj.MinActivityAtSpecificPx = min(ThisObj.DataMatrix(yk,xk,:));
            catch
                ThisObj = ThisObj.GetVarMap;
                [row,col] = find(ThisObj.VarMap == max(max(ThisObj.VarMap)));
                ThisObj.MaxActivityAtSpecificPx = max(ThisObj.DataMatrix(row,col,:));
                ThisObj.MinActivityAtSpecificPx = min(ThisObj.DataMatrix(row,col,:));
            end
            
            % Initialize Monolayer structure
            ThisObj.Monolayer.radius_cm = 1;
            ThisObj.Monolayer.radius_px = [];
        end
                
        function SensorA = GetIBI(SensorA, SpaceTime)
            % Not very accurate. Fix ComputeIBI_IndividualThreshold().
            SensorA = SensorA.GetFullTimeSeries(SpaceTime);
            SensorA = SensorA.GetAPCrossTimes(SpaceTime);
            SensorA.IBI = ComputeIBI_IndividualThreshold(SensorA, SpaceTime);
        end
        
        function SensorA = GetAPCrossTimes(SensorA, SpaceTime)
            [SensorA.APCrossTimes, SensorA.APThreshold] = GetCrossTimes(SensorA, SpaceTime);
        end
        
        function SensorA = GetTimeSeries(SensorA, SpaceTime)
            % First of all, clear TimeSeries to store new data
            SensorA.TimeSeries = [];
            % Extract interval if FullTimeSeries is already available
            if isempty(SensorA.FullTimeSeries) == 0
                
                taxis = SensorA.TimeAxis;      % Time axis in sec
                
                if SpaceTime.ti > SensorA.TotalFrames*SensorA.SampPeriod || ... 
                        SpaceTime.tf > SensorA.TotalFrames*SensorA.SampPeriod
                    error('Time outside bounds.')
                else
                    ind = taxis >= SpaceTime.ti & taxis <= SpaceTime.tf;
                    SensorA.TimeSeries(:,1) = taxis(ind);
                    for i = 1 : size(SensorA.FullTimeSeries,2)
                        SensorA.TimeSeries(:,i+1) = SensorA.FullTimeSeries(ind, i);
                    end
                end
            else 
                % ...if not, compute full time series and call function again
                SensorA = SensorA.GetFullTimeSeries(SpaceTime);
                SensorA = SensorA.GetTimeSeries(SpaceTime);
            end
            
        end
        
        function SensorA = GetFullTimeSeries(SensorA, SpaceTime)
            SensorA.FullTimeSeries = ComputeTimeSeries(SensorA, SpaceTime);
        end
        
        function SensorA = GetxCorr(SensorA, SpaceTime)
            [SensorA.xCorr, SensorA.lag] = ComputePartCorr(SensorA, SpaceTime);
        end
        
        function SensorA = GetVarMap(SensorA,varargin)
            % Computes the variance map of movie using all frames.
            if isempty(SensorA.VarMap)==1
                [~,SensorA.VarMap] = svmMAP(SensorA.DataMatrix,varargin{:}); 
            end         
        end
        
        function SensorA = GetRadius(SensorA, varargin)
            % If SensorA.Monolayer.radius_px is empty, computes monolayer
            % radius from the two points selected by user in the variance
            % map that is plotted; else, ignores request. To re-do the
            % operation, use the flag 'redo' as argument.
            % Example:
            % mysensor.GetRadius('redo');
            %
            flag_redo = ismember('redo', varargin);
            if isempty(SensorA.Monolayer.radius_px)==1 || flag_redo==1
                SensorA = SensorA.GetVarMap;
                pcolor(SensorA.VarMap); shading interp
                axis square
                disp('Select initial and final points on image.')
                [x,y] = ginput(2);
                SensorA.Monolayer.radius_px = sqrt( (x(2)-x(1))^2 + (y(2)-y(1))^2 )/2; 
                SensorA.Monolayer.radius_px = round(SensorA.Monolayer.radius_px);
            end
        end
        
        %% getter function for sampling period
        function SampPeriod = get.SampPeriod(SensorA)
            SampPeriod = 1/SensorA.SampFreq;
        end
        
        %% getter function total number of frames
        function TotalFrames = get.TotalFrames(SensorA)
            TotalFrames = size(SensorA.DataMatrix,3);
        end
        
         %% getter function for total time (in s)
        function TotalTime = get.TotalTime(SensorA)
            TotalTime = SensorA.TotalFrames * SensorA.SampPeriod;
        end
        
         %% getter function for time axis (in s)
        function TimeAxis = get.TimeAxis(SensorA)
            TimeAxis = (1 : SensorA.TotalFrames) * SensorA.SampPeriod;
        end
        
         %% getter function for FLAG_STIM
        function FLAG_STIM = get.FLAG_STIM(SensorA)
            FLAG_STIM = CheckStim(SensorA);
        end
        %% 
        % Action Potential crossing times and thresholds for each AP
        function [crosstimes, y_th] = GetCrossTimes(SensorA, SpaceTime)
            for i = 1 : SpaceTime.nCoordinates
                [crosstimes{i}, y_th{i}] = ComputeCrossTimes(SensorA.TimeAxis, SensorA.FullTimeSeries(:,i));
            end       
        end
        
        % Stimuli crossing times
        function StimRef = ComputeStimReference(SensorA)
            if SensorA.FLAG_STIM==1
                t = SensorA.TimeAxis;
                y = SensorA.Stimuli;
                y_th = SensorA.StimuliThreshold;
                sign_th = 1;
                StimRef = CrossTimes_SharedThreshold(t, y, y_th, sign_th);
            else
                StimRef = [];
            end
        end
        
        % Get the time series for the full interval.
        % For Type '2-coord':
        % First column : time series of the initial pixel (xi, yi). 
        % Second column: time series of the final   pixel (xf, yf).
        % For Type 'n-coord':
        % k-th column : time series of the k-th pixel (xk, yk). 
        function TimeSeries = ComputeTimeSeries(SensorA, SpaceTime)
            if strcmp(SpaceTime.SetupType, '2-coord')==1
                TimeSeries(:,1) = SensorA.DataMatrix(SpaceTime.yi, SpaceTime.xi,:);
                TimeSeries(:,2) = SensorA.DataMatrix(SpaceTime.yf, SpaceTime.xf,:);
            else
                TimeSeries = zeros([SensorA.TotalFrames, SpaceTime.nCoordinates]);
                j = 1;
                for q = 1 : SpaceTime.nCoordinates
                    TimeSeries(:,j) = SensorA.DataMatrix(SpaceTime.xy(2*q), SpaceTime.xy(2*q-1),:);
                    j = j + 1;
                end
            end
%             if strcmp(SpaceTime.SetupType, '2-coord')==1
%                 TimeSeries(:,1) = AvgPx(SensorA, SpaceTime.xi, SpaceTime.yi);
%                 TimeSeries(:,2) = AvgPx(SensorA, SpaceTime.xf, SpaceTime.yf);
%             else
%                 TimeSeries = zeros([SensorA.TotalFrames, SpaceTime.nCoordinates]);
%                 j = 1;
%                 for q = 1 : SpaceTime.nCoordinates
%                     TimeSeries(:,j) = AvgPx(SensorA, SpaceTime.xy(2*q-1), SpaceTime.xy(2*q));
%                     j = j + 1;
%                 end
%             end
        end
        
        % UPDATE THIS FUNCTION
        % Get the time series for the interval specified in SpaceTime.
        % First column : time series of the initial pixel (xi, yi). 
        % Second column: time series of the final   pixel (xf, yf).
        % Third column : time axis.
        function TimeSeries = ExtractTimeSeries(SensorA, SpaceTime)
            if strcmp(SpaceTime.SetupType, '2-coord')==1
                [TimeSeries(:,1), ...
                 TimeSeries(:,2), ...
                 TimeSeries(:,3)] = ExtractInterval(SensorA, SpaceTime);
            end
        end
        
        %
        function IBI = GetIBI_SharedThreshold(SensorA, SpaceTime)
           IBI = ComputeIBI_SharedThreshold(SensorA, SpaceTime);
        end
        
        % Compute the wave velocity for the time interval specified in
        % SpaceTime
        function WaveVel = ComputeWaveVel(SensorA, SpaceTime)
            if strcmp(SpaceTime.SetupType, '2-coord')==1
                WaveVel = WaveVelocity(SensorA, SpaceTime);
            else
                error('SpaceTime setup not compatible. Use 2-coord setup.')
            end
        end
        
        % Compute the cross-correlation for positive lags between the two
        % time series for the full interval
        function [xCorr, lag] = ComputeFullCorr(SensorA, SpaceTime)
            Y = ComputeTimeSeries(SensorA, SpaceTime);
            [xCorr, lag] = GetCorr(Y(:,2),Y(:,1));
            lag = lag' * SensorA.SampPeriod;
        end
        
        % Compute the cross-correlation for positive lags between the two
        % time series within the time interval specified in SpaceTime
        function [xCorr, lag] = ComputePartCorr(SensorA, SpaceTime)
            Y = ExtractTimeSeries(SensorA, SpaceTime);
            [xCorr, lag] = GetCorr(Y(:,2),Y(:,1));
            lag = lag' * SensorA.SampPeriod;
        end
        
        % Plot graphs related to the method used to get velocity
        function plotVMethod(SensorA)
            figure
            subplot 211
            plot(SensorA.TimeSeries(:,3), SensorA.TimeSeries(:,1), ...
                 SensorA.TimeSeries(:,3), SensorA.TimeSeries(:,2))
            legend('Initial px', 'Final px')
            xlabel 'Time (s)'
            ylabel 'Fluorescence (AU)'
            
            subplot 212
            plot(SensorA.lag, SensorA.xCorr)
            xlabel 'Time (s)'
            ylabel 'Cross-Correlation'
        end
    end
    
end

%% Auxiliary functions to this class
% n-coord type. Individual threshold ==================================== %
function IBI = ComputeIBI_IndividualThreshold(SensorA, SpaceTime)
% Stores the IBI values and the IBI_times of each coordinate in columns of
% each cell. Column 1 is time (s), column 2 is IBI (s).

    ncoord = length(SpaceTime.xy)/2;
    IBI = cell(1,ncoord);
    for i = 1 : ncoord
        if isempty(SensorA.APCrossTimes) == 0
            tcross_up = SensorA.APCrossTimes{i};
        else
            % Saves the time when crosses threshold
            [tcross_up, ~] = ComputeCrossTimes(SensorA.TimeAxis, SensorA.FullTimeSeries(:,i));  
        end
        ibi = diff(tcross_up);                            
        IBI{i}(:,2) = ibi;
        IBI{i}(:,1) = tcross_up(1:end-1) + ibi/2;
    end
end

function [crosstimes, y_th] = ComputeCrossTimes(t, y)

    % initialize height of individual AP
    delta_y = 0;        
    
     % height of AP to consider that a peak is an AP
    [~, delta_y_th] = GetAPHeights(y);
    
    % percentage of height of an AP to consider as crossing threshold
    PercOfAP = 0.5;
    
    FLAG_VALE = 0;
    FLAG_PEAK = 0;

    j = 1;              % iterator for the upstroke array
    it = 1;             % iterator for the crosstimes and y_th arrays:
                        % the upstroke array resets at every valley, and
                        % each AP has a threshold y_th.
                        % upstroke(:,1) -> values of y for an AP
                        % upstroke(:,2) -> values of t for the same AP
    
    % n points can form at most (n-1)/2 peaks                    
    npeaks = floor( (length(y)-1)/2 );
    nPointsPerUpstroke = 100;
    
    crosstimes = zeros([npeaks, 1]);
    y_th       = zeros([npeaks, 1]);  
    upstroke   = zeros([nPointsPerUpstroke,2]);
    
    for q = 2 : length(y)-1
        if ( y(q-1) >= y(q) && y(q+1) > y(q) ) % at valley, reset upstroke and j
            FLAG_VALE = 1;
            FLAG_PEAK = 0;
            upstroke = [];
            j = 1;
        end

        if y(q-1) < y(q) && y(q+1) < y(q) && FLAG_VALE == 1 % at peak after valley
            upstroke(j,1) = y(q);
            upstroke(j,2) = t(q);
            j=j+1;
            FLAG_PEAK = 1;
            FLAG_VALE = 0;
         elseif y(q-1) < y(q) && y(q) < y(q+1) && FLAG_PEAK == 0 && FLAG_VALE == 1 % between valley and peak         
            
            if rad2deg( atan2(y(q)-y(q-1),t(q)-t(q-1)) ) > 70
                upstroke(j,1) = y(q);
                upstroke(j,2) = t(q);
                j=j+1;
            end

        end
        
        if FLAG_PEAK==1
            delta_y = abs(upstroke(end,1) - upstroke(1,1));
            FLAG_PEAK = 0;
        end
        
        if delta_y >= delta_y_th
%             plot(upstroke(:,2),upstroke(:,1), 'o'); hold on;
            for k = 2 : length(upstroke)                
                y_th(it) = upstroke(1,1) + PercOfAP * delta_y;  % threshold  
                
                if upstroke(k-1,1) < y_th(it) && upstroke(k,1) >= y_th(it) % if crosses threshold
                    crosstimes(it) = interp1(upstroke(:,1), upstroke(:,2), y_th(it));
                    it = it + 1;
                end   
            end
            delta_y = 0;    % Resets delta_y
        end
    end
    
    % Removes the excedent allocated space
    crosstimes(crosstimes==0) = [];
    y_th(y_th == 0) = [];
    
    % Removes the repetead element (due to where 'it' was incremented)
    y_th(end) = [];
end

function [APHeights, APHeights_th] = GetAPHeights(y)
% [APHeights, APHeights_th] = GetAPHeights(y) computes the height of all
% peaks in a time series, including noise. Then it builds an histogram of
% the height values and uses kmeans to cluster this distribution into 3
% groups. The x-value of the central group is an estimate of the threshold
% for action potential.

it = 1; % iterator for delta_y

% n points can form at most (n-1)/2 peaks                    
npeaks = floor( (length(y)-1)/2 );  
delta_y    = zeros([npeaks, 1]);   

% estimating 2^7 points max for an upstroke
upstroke   = zeros([2^7,2]);        

% Flags to identify peaks and valleys
FLAG_VALE = 0;
FLAG_PEAK = 0;

j = 1;      % iterator for the upstroke array 

for q = 2 : length(y)-1
        if ( y(q-1) >= y(q) && y(q+1) > y(q) ) % at valley, reset upstroke and j
            FLAG_VALE = 1;
            FLAG_PEAK = 0;
            upstroke = [];
            j = 1;                         
        end

        if ~( y(q-1) < y(q) && y(q+1) < y(q) ) && FLAG_VALE == 1 % between valley and peak         
            upstroke(j) = y(q);
            j=j+1;
        elseif y(q-1) < y(q) && y(q+1) < y(q) && FLAG_VALE == 1  % at peak, after a valley
            upstroke(j) = y(q);
            j=j+1;
            FLAG_PEAK = 1;
            FLAG_VALE = 0;
        end
        
        if FLAG_PEAK==1 % at peak, compute height
            delta_y(it) = upstroke(end) - upstroke(1);
            FLAG_PEAK = 0;
            it = it + 1;
        end
end
delta_y(delta_y==0) = [];

[counts, centers] = hist(delta_y); 

X = [centers' counts'];

% Clustering the values of delta into 3 groups: spikes, noise and
% 'transition'. The ctrs of the clusters are the mean of the heights. So
% we're interested in the maximum value of ctrs as the AP height, and in
% the middle value of ctrs, which should be a good estimate of the height
% threshold.
k = [2 3 4];
s_avg = zeros([length(k),1]);
for i = 1 : length(k)
    [idx,ctrs] = kmeans(X,k(i),...
                'Distance','sqeuclidean',...
                'Replicates',5);
    s = silhouette(X, idx, 'sqeuclidean');
    s_avg(i) = mean(s);
end
best_k = k( s_avg==max(s_avg) );

switch best_k
    case 2
        % Get th as mean
        APHeights_th = (ctrs(1,1) + ctrs(2,1))/2;
    case 3
        % Get the threshold as the intermediate value
        for i = 1 : 3
            if ctrs(i,1) > 1.1*min(ctrs(:,1)) && ctrs(i,1) < 0.9*max(ctrs(:,1))
                APHeights_th = ctrs(i,1);
            elseif abs(min(ctrs(:,1)) - max(ctrs(:,1))) < 0.5
                APHeights_th = ctrs(1,1)/2;
            end
        end
    case 4
        SortCtrs = sort(ctrs(:,1));
        APHeights_th = (SortCtrs(2)+SortCtrs(3))/2;
end              
APHeights = max(ctrs(:,1));

%------------------------------------------------------------------------ %
% Uncoment the lines below to check the distribution of height values
figure 
subplot 221
plot(delta_y,'o');
title 'Height values'
ylabel 'Height'
xlabel 'Point #'

subplot 222
hist(delta_y)
title 'Height distribution'
ylabel 'Count'
xlabel 'Height'

subplot 223
[silh1, ~] = silhouette(X, idx, 'sqeuclidean');
% Vertical line at the mean value
hl = line([mean(silh1) mean(silh1)], ylim);
hl.Color = 'red';
hl.LineWidth = 1.5;
title 'Silhouette plot'

subplot 224
hp1 = plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12);
hold on
hp2 = plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12);
hp3 = plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',12);
plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
       'Location','NE');
% ----------------------------------------------------------------------- %

if APHeights_th < 0.20 * APHeights
    error('APHeights_th < 0.2APHeight. Double check clustering in file Sensor.m, function GetAPHeights(y).')
else
    fprintf('k = %d \t AP mean height: %.2f   //   '  , best_k, APHeights);
    fprintf('AP height threshold: %.2f \n', APHeights_th);
end

end
%% Shared threshold ===================================================== %
function IBI = ComputeIBI_SharedThreshold(SensorA, SpaceTime)
% Stores the IBI values and the IBI_times of each coordinate in columns of
% each cell. Column 1 is time (s), column 2 is IBI values (s).
    
    % Ideally, eps == refractory period duration
    eps = 0.20; % (in seconds)
    ncoord = length(SpaceTime.xy)/2;
    IBI = cell(1,ncoord);
    for i = 1 : ncoord
        % Saves the time when crosses threshold
        tcross_up_raw = CrossTimes_SharedThreshold(SensorA.TimeAxis, SensorA.FullTimeSeries(:,i), SensorA.APThreshold, 1);    
        tcross_up     = CleanCrossTimes(tcross_up_raw, eps);
        ibi = diff(tcross_up);
        IBI{i}(:,2) = ibi;
        IBI{i}(:,1) = tcross_up(1:end-1) + ibi/2;
    end
    
end

function [tcross, ind] = CrossTimes_SharedThreshold(t, y, y_th, sign_th)
    % Min-Young Kim, June 8, 2006
    % To detect crossing points from below to above threshold in real vector y
    % event triggering: y = y_th, with sign( y'(y_th) ) = sign_th 
    % sign_th = 1 for positive slope, -1 for negative slope
    % t_crossing is estimated using linear interpolation
    % If (nx) is one of the crossing point, 
    % t_crossing(nx) = t(nx-1) + (y_th-y(nx-1)/(y(nx) - y(nx-1))*delta_t
    
    [~, n_coly] = size(y);
    [~, n_colt] = size(t);
    if n_coly == 1
        y = y';
    end            % making a row vector
    
    if n_colt == 1
        t = t';
    end           
    
    dy = sign_th*(y-y_th);
    ind = find( ([0 dy]<0) & ([dy 0]>=0)); % indice right after crossing 
    ind = ind(1:end-1);     % remove [the first and] the last crossing

    delta_t= t(2)-t(1);
    tcross = t(ind-1) + (y_th - y(ind-1))./(y(ind) - y(ind-1))*delta_t;    
end

function tcross = CleanCrossTimes(tcross, eps)
% function tcross = clean_tcross(tcross, eps) cleans the array tcross so
% that the difference tcross(k)-tross(k-1) is always greater than eps.
    IBI = diff(tcross);
    temp = isempty(find(IBI <= eps, 1));

    if temp == 0
        for i = 2 : numel(tcross)
            tempdiff = tcross(i) - tcross(i-1);
            if tempdiff <= eps
                tcross(i) = [];
                break
            end
        end
    else
        return
    end
    tcross = clean_tcross(tcross, eps);
end

%% 2-coord type. 'Initial and final'-related ============================ %
function [xCorr, lag] = GetCorr(a, b)    
    [xCorr, lag] = xcorr(a, b);
%     ind = find(lag >=0);
%     xCorr = xCorr(ind);
%     lag = lag(ind);
end

function [Yi_partial, Yf_partial, taxis] = ExtractInterval(SensorA, SpaceTime)
% Returns the time series of the initial and final pixel and the time axis
% within the specified interval
    Yi = AvgPx(SensorA, SpaceTime.xi, SpaceTime.yi);
    Yf = AvgPx(SensorA, SpaceTime.xf, SpaceTime.yf);

    % Time axis in sec
    taxis = SensorA.TimeAxis; % sec

    % Extract interval
    if SpaceTime.ti > SensorA.TotalFrames*SensorA.SampPeriod || ... 
            SpaceTime.tf > SensorA.TotalFrames*SensorA.SampPeriod
        error('Time outside bounds.')
    else
        ind = find(taxis >= SpaceTime.ti & taxis <= SpaceTime.tf);
        Yi_partial = Yi(ind);
        Yf_partial = Yf(ind);
        taxis = taxis(ind);
    end
end

function wave_vel = WaveVelocity(SensorA, SpaceTime)
    % Extract interval
    [Yi_partial, Yf_partial, ~] = ExtractInterval(SensorA, SpaceTime);

    % Cross correlation
    [xCorr_partial, lag_partial] = GetCorr(Yf_partial, Yi_partial);
    lag_partial = lag_partial * SensorA.SampPeriod;

    % where C is max
    indmax = xCorr_partial == max(xCorr_partial);
    wave_vel = SpaceTime.Distance/lag_partial(indmax);
end

%% Stimuli-related ====================================================== %
function FLAG_STIM = CheckStim(SensorA)                 % Boolean
% FLAG_STIM = CheckStim(SensorA) checks if the signal from "Stim" is
% from actual stimuli recording or not.

    t = SensorA.TimeAxis;                    % It's already in seconds
    threshold = SensorA.StimuliThreshold;                        % max(stim) ~ 10^4

    tcross_up = CrossTimes_SharedThreshold(t, SensorA.Stimuli, threshold, 1);

    d = diff(tcross_up);                    % ISI
    if numel(d) < 3
        FLAG_STIM = 0;
    else
        FLAG_STIM = 1;
    end
end 

function [StimPeriod, StimRange] = GetStimInfo(SensorA) % (in seconds)
% Function to get the intervals in which the stimuli happen, and the period
% of each stimulation. Each line of StimRange refers to a group; the first
% column is the time when that group started and the second column is the
% time when it ended. 

    k = 2;                                   % Percentage of 1st isi to divide into groups
    t = SensorA.TimeAxis;                    % It's already in seconds
    threshold = SensorA.StimuliThreshold;    % To identify spikes (10% of the max)
    tcross_up = CrossTimes_SharedThreshold(t, SensorA.Stimuli, threshold, 1);
    
% --- Another way of getting tcross ------------------------------------- %
%     [~, ind_ac] = CrossTimes_SharedThreshold(t, SensorA.Stimuli, threshold, 1);
%     y_b = SensorA.Stimuli(ind_ac-1);
%     y_a = SensorA.Stimuli(ind_ac);
%     t_b = t(ind_ac-1);
%     t_a = t(ind_ac);
%     y = [y_b; y_a];
%     tt= [t_b; t_a];
%     for i = 1 : length(y_b);
%         tcross_up(i) = interp1(y(:,i),tt(:,i),threshold);
%     end
% ----------------------------------------------------------------------- %

    d = diff(tcross_up);                    % ISI
    
    if SensorA.FLAG_STIM==0
        StimPeriod = [];
        StimRange  = [];
        return;
    else
        index = find( d > k*d(1) );         % when ISI > first ISI value
    end
  
    ndiv = length(index);                   % Number of dividers
    ngroups = ndiv + 1;                     % Number of groups
    
    StimPeriod = zeros([ngroups, 1]);       % Stores the mean values
    StimRange  = zeros([ngroups, 2]);       % Stores the start (1) and end (2) time for each group
    Groups = cell(1,ngroups);               % Stores the values of 'd' for each group
    ind = 1;                                % Index to define the Groups range
    
    % Indexes for the StimRange
    ind1 = 1;
    if ngroups==1
        ind2 = length(tcross_up);
    else
        ind2 = index(1);
    end
    
    for i = 1 : ngroups
        if ngroups-i == 0 % Last group
            Groups{i} = d(ind:length(d));
            StimRange(i,1) = tcross_up(ind1);
            StimRange(i,2) = tcross_up(ind2);        
        else
            Groups{i} = d(ind:index(i)-1);
            StimRange(i,1) = tcross_up(ind1);
            StimRange(i,2) = tcross_up(ind2);
            ind = index(i)+1;

            if i+1 <= length(index) % Same as 'if i <= ngroup'
                ind1 = ind2 + 1;
                ind2 = index(i+1);
            else
                ind1 = ind2 + 1;
                ind2 = length(d)+1;            
            end
        end
        StimPeriod(i,1) = round(mean(Groups{i}), 3); % round to 3 decimal places
    end
end
