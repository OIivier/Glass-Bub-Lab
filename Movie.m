% Class used to generate movie from an object of type Sensor. 
%
% USAGE: Movie(Some_Sensor_Data, optional_arguments). 
% If an optional argument is not used, the default value will be used.
% Below is a list of the optional arguments and their default values.
%
% OPTIONAL ARGUMENTS: 
%       initialframe    - first frame of movie               default: 1
%       finalframe      - last frame of movie                default: max frame of data
%       speed           - screening speed                    default: 1
%       frmDelay        - # of periods to delay/frame        default: 1
%       save            - true: save movie as gif            default: false         
%       pause           - true: press buttom to pass movie   default: false          
%       xy              - [x1,y1,...,xn,yn] plots (xi,yi)    default: [] 
%
% Examples:
% 1.      
%       Data = Sensor(SF1, Stim);
%       Movie(Data, ...
%            'initialframe', 1000, ...
%            'finalframe'  , 2000, ...
%            'speed'       , 2,    ...
%            'xy'          , [50, 20, 25, 40])
% 
% This code will generate a movie from frame 1000 to frame 2000 with
% screening speed 2 and plot the coordinates (50,20) and (25,40).
%
% 2.      
%       Data = Sensor(SF1, Stim);
%       Movie(Data)
% 
% This code will generate a movie using the default value for all
% arguments: from frame 1 to the last frame of SF1, with screening speed 1.
%
%
% L. Campanari, Feb 2016.

classdef Movie
    
    properties
        InitialFrame
        FinalFrame
        Speed                           % Screening speed
        FrmDelay                        % # of periods (0.025 ms) the gif should delay at each frame
        Save
        Pause
        XY                              % Coordinates to plot on the movie
        hfig
        ax
        Title
        Position
        FileName
        Xlabel
        Ylabel
        Axis
        XTick
        YTick
        Colormap
    end
    
    methods
        
        function ThisMovie = Movie(SensorA, varargin)
           p = inputParser;          
           % ------------------------------------------------------------ %
           % Default values and validation
           defaultInitialFrame = 1;
           checkInitialFrame = @(x) isnumeric(x) && (x >= 1) && (x < SensorA.TotalFrames);
           
           defaultFinalFrame = SensorA.TotalFrames;
           checkFinalFrame = @(x) isnumeric(x) && (x > 1) && (x <= SensorA.TotalFrames);
           
           defaultSpeed = 1;
           checkSpeed = @(x) isnumeric(x) && (x > 0);
           
           defaultFrmDelay = 1;
           checkFrmDelay = @(x) isnumeric(x) && (x >= 0);
           
           defaultSave = 'false';
           validSave = {'true', 'false'};
           checkSave =  @(x) any(validatestring(x, validSave));
           
           defaultPause = 'false';
           validPause = {'true', 'false'};
           checkPause =  @(x) any(validatestring(x, validPause));
           
           defaultTitle = '';
           checkTitle = @(x) 1==1;
           
           defaultPos = [];
           checkPos = @(x) isnumeric(x) && numel(x)==4;
           
           defaultFileName = [];
           checkFileName = @(x) 1==1;
           
           defaultXlabel = 'Pixel';
           checkXlabel = @(x) 1==1;
           
           defaultYlabel = 'Pixel';
           checkYlabel = @(x) 1==1;
           
           defaultAxis = 'on';
           validAxis = {'on', 'off'};
           checkAxis =  @(x) any(validatestring(x, validAxis));
           
           defaultColormap = 'summer';
           checkColormap = @(x) 1==1;
              
           % Add previously defined stuff to inputParser
           p.addRequired('SensorA', @(x) isa(x, 'Sensor'));
           p.addParameter('initialframe', defaultInitialFrame, checkInitialFrame   ); % 1
           p.addParameter('finalframe',   defaultFinalFrame,   checkFinalFrame     ); % 2
           p.addParameter('speed',        defaultSpeed,        checkSpeed          ); % 3
           p.addParameter('xy',           [],              @(x) mod(length(x),2)==0); % 4
           p.addOptional('pause',         defaultPause,        checkPause          ); % 5
           p.addParameter('save',         defaultSave,         checkSave           ); % 6 
           p.addParameter('frmdelay',     defaultFrmDelay,     checkFrmDelay       ); % 7
           p.addParameter('title',       defaultTitle,         checkTitle          ); % 8
           p.addParameter('pos',         defaultPos,           checkPos            ); % 9
           p.addParameter('filename',    defaultFileName,      checkFileName       ); % 10
           p.addParameter('xlabel',      defaultXlabel,        checkXlabel         ); % 11
           p.addParameter('ylabel',      defaultYlabel,        checkYlabel         ); % 12
           p.addParameter('axis',        defaultAxis,          checkAxis           ); % 13
%            p.addParameter('XTick',        defaultXTick,        checkXTick          ); % 14
%            p.addParameter('YTick',        defaultYTick,        checkYTick          ); % 15
           p.addParameter('colormap',    defaultColormap,      checkColormap           ); % 13 
           % ------------------------------------------------------------ %
           
           parse(p, SensorA, varargin{:}); 
           
           % Attributes parameters to class properties
           ThisMovie.InitialFrame = p.Results.initialframe;
           ThisMovie.FinalFrame   = p.Results.finalframe;
           ThisMovie.Speed        = p.Results.speed;
           ThisMovie.Save         = p.Results.save;
           ThisMovie.XY           = p.Results.xy;
           ThisMovie.Pause        = p.Results.pause;
           ThisMovie.FrmDelay     = p.Results.frmdelay;
           ThisMovie.Title        = p.Results.title;
           ThisMovie.Position     = p.Results.pos;
           ThisMovie.FileName     = p.Results.filename;
           ThisMovie.Xlabel       = p.Results.xlabel;
           ThisMovie.Ylabel       = p.Results.ylabel;
           ThisMovie.Axis         = p.Results.axis;
%            ThisMovie.XTick        = p.Results.XTick;
%            ThisMovie.YTick        = p.Results.YTick;
            ThisMovie.Colormap    = p.Results.colormap;
           
           % Call movie maker
           ThisMovie = ThisMovie.Producer(SensorA);
        end
        
        function MyMovie = Producer(MyMovie, Sensor)
            [MyMovie.hfig, MyMovie.ax] = MakeMovie(MyMovie, Sensor);
        end
    end
end

%% Auxiliary functions to this class
function [hfig, ax] = MakeMovie(MyMovie, Sensor)
tic
close all;

if isempty(Sensor.Stimuli)==1
    FLAG_STIM = 0;
else
    % To plot frequency and stimulus signal
    FLAG_STIM = 1;
    groupk = 1;
    threshold_up = 1000; %Threshold to consider it's a stimulus
end

% Sampling period
T  = Sensor.SampPeriod;        

PositiveSave = {'save', 'true'};
ValidateSave = any(strcmp(MyMovie.Save, PositiveSave));

PositivePause = {'pause','true'};
ValidatePause = any(strcmp(MyMovie.Pause, PositivePause));

% Array of frame numbers
frame = MyMovie.InitialFrame : MyMovie.Speed : MyMovie.FinalFrame;
nFrm = length(frame);

% Allocate structure F for getframe();
F = struct('cdata', cell(1,nFrm), 'colormap', cell(1,nFrm)); 
hfig = figure;

try
    colormap(MyMovie.Colormap)
catch
    colormap(summer)
    disp('Colormap not reconized. Using ''summer''.')
end

if isempty(MyMovie.Title)
    if ~isempty(Sensor.Name)
        titlestr = ['File: ' Sensor.Name '     ' 'Date: ' Sensor.Date];
    else
        titlestr = '';
    end
else
    titlestr = MyMovie.Title;
end
    
for i=1:length(frame)-1
    
    ax = pcolor(Sensor.DataMatrix(:,:,frame(i)));
    shading interp % original: interp
    cmin=Sensor.MinActivityAtSpecificPx; 
    cmax=Sensor.MaxActivityAtSpecificPx;    
    caxis([cmin cmax]);
    title(titlestr);
    hold on;
    if FLAG_STIM ~= 0
        PlotStim(Sensor.Stimuli(frame(i)), Sensor.Stimuli(frame(i+1)) ,threshold_up);
    end
    
    axis square;
    axis(MyMovie.Axis)
    xlabel(MyMovie.Xlabel)
    ylabel(MyMovie.Ylabel)
    
    set(gcf,'color','w');
    set(gca,'YDir','reverse');
    set(gca,'layer','top', ...
            'linewidth',1,...
            'FontSize',11,...
            'FontWeight','normal');
    
    if ~isempty(MyMovie.Position)
        set(hfig,'pos',MyMovie.Position)    
    end
    
    % Plot specified points
    for q = 1 : 2 : length(MyMovie.XY)
        plot(MyMovie.XY(q), MyMovie.XY(q+1), 'w*')
        text(MyMovie.XY(q), MyMovie.XY(q+1)-3, ...
            sprintf('%d',(q+1)/2), ...  
            'Color', 'w',                 ...
            'FontSize', 11,               ...
            'HorizontalAlignment','center',...   
            'VerticalAlignment','top');   
    end
        
    % Print time in seconds
    text(77,2,[sprintf('%.2f',frame(i)*T) 's'],      ...   % <--
        'Color', 'w',                 ...
        'FontSize', 11,               ...
        'FontWeight', 'normal',         ...
        'HorizontalAlignment','right',...   
        'VerticalAlignment','cap');
    
    % Print current value of frequency of stimulation
    if FLAG_STIM ~= 0
        if frame(i)*T >= Sensor.StimuliRange(groupk,1) && ...
            frame(i)*T <= Sensor.StimuliRange(groupk,2)
               text(77,5,sprintf( '%d ms', ...
                   10*round(100*Sensor.StimuliPeriod(groupk)) ), ...
                    'Color', 'w',                 ...
                    'FontSize', 11,               ...
                    'HorizontalAlignment','right',...   
                    'VerticalAlignment','cap');
        elseif frame(i)*T > Sensor.StimuliRange(groupk,2)
            if groupk < size(Sensor.StimuliRange,1)
                groupk = groupk+1;
            end
        end
    end
    hold off;
    F(i)=getframe(gcf);
    
    % Wait for buttom press?
    if ValidatePause ~= 0
        waitforbuttonpress                    
    end
    
    % Save movie?
    if ValidateSave ~= 0
        if ~isempty(MyMovie.FileName)
            gifname = MyMovie.FileName;
        else
            if isempty(Sensor.Name) == 0
                if strcmp(Sensor.Name(end-2:end),'.da')
                    gifname = [Sensor.Name(1:end-3) '_' Sensor.Date];
                else
                    gifname = [Sensor.Name '_' Sensor.Date];
                end
            else
                gifname = 'mygif';
            end
        end
        
        drawnow
        im = frame2im(F(i));
        [A,map] = rgb2ind(im,256); 
        if i == 1;
            imwrite(A,map,[gifname '.gif'],'gif','LoopCount',Inf,'DelayTime',T*MyMovie.FrmDelay);
        else
            imwrite(A,map,[gifname '.gif'],'gif','WriteMode','append','DelayTime',T*MyMovie.FrmDelay);
        end
    end
end

toc
end

function PlotStim(Stim0, Stim1, th)

threshold_up = th;
threshold_down = -threshold_up;
posx = 62;
posy = 7;
size = 54;

if Stim0 < threshold_up && Stim1 >= threshold_up
    plot(posx,posy, 'g.', 'markers', size)
end

% if Stim0 > threshold_down && Stim1 <= threshold_down
%     plot(posx,posy, 'y.', 'markers', size)
% end
    
end


