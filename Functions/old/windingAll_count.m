% Last modified on 12.04.17

% This algorithm detects the position of phase singularities based on the 
% calculation of the winding number inside a given circle of given radius.

function [NewSourceMap] = windingAll_count(phaseMap, iniFrame, dura, speed, save)
% Maps the points that give a winding number of -1 or +1 on the frame. 
% White = clockwise spiral (winding# = -1)
% Black = counterclockwise spiral (winding# = +1)

tic
%% Initializing variables
temp = phaseMap+180; % shifts all values up 180 so that the angles run from 0 to 360 rather than -180 to 180.  
[row, column, frame] = size(temp);
video = phaseMap;
sourceMap = zeros(row, column, 2);

% This varies the radius depending on the length and height of phasemap. 
radius = round(row/10);
if radius == 0
    radius = 1;
end 

radius_check = radius
%% Variables for the video 
Fs = 40;              % Sampling frequency
T = 1/Fs;             % Sampling period

finFrame = iniFrame + dura; % sets the duration to include a certain amount of cycles 
% finFrame = ceil(iniFrame + (Period*Fs*10)); % sets the duration to include a certain amount of cycles 
cycle = 1;

finFrame
%% WINDING NUMBER

for t = iniFrame:speed:finFrame
% For saving when winding number = -1 or +1 
clocktic = 1;
countertic = 1; 
clock = zeros(1,1);
counterClock = zeros(1,1);
%% Draw a circle around the center and calculate the winding number 
for centerX = 1:column
    
    % Make sure that the circle doesn't go out of the frame. 
    if (centerX - radius)<1 || (centerX + radius)>column 
        continue % continues the loop  
    end 

    for centerY = 1:row
        
        % Make sure that the circle doesn't go out of the frame. 
        if(centerY - radius)<1 || (centerY + radius)>row 
            continue 
        end 

    % Reset these values for every center point. 
    count = 1;
    arrayCirc = zeros(1,1,frame); % arrayCirc saves the values on the circle at ALL FRAMES  
    windingNum = 0;

    % Draw a circle around (centerX, centerY) and save the values in the circle to arrayCirc 
    % Quadrant 1
    for x = centerX:column  
        for y = 1:centerY
            r = round(sqrt((centerY-y).^2+(centerX-x).^2));
            if r == radius
               arrayCirc(1,count,:) = temp(y,x,:);
               count = count + 1; 
            end
        end 
    end

    % Quadrant 4 
    for x = column:-1:centerX
        for y = centerY+1:row
            r = round(sqrt((centerY-y).^2+(centerX-x).^2));        
            if r == radius 
               arrayCirc(1,count,:) = temp(y,x,:);
               count = count + 1; 
            end
        end
    end 

    % Quadrant 3  
    for x = (centerX-1):-1:1  
        for y = row:-1:centerY+1
            r = round(sqrt((centerY-y).^2+(centerX-x).^2));        
            if r == radius 
               arrayCirc(1,count,:) = temp(y,x,:);
               count = count + 1; 
            end
        end
    end

    % Quadrant 2 
    for x = 1:centerX-1  
        for y = centerY:-1:1
            r = round(sqrt((centerY-y).^2+(centerX-x).^2));        
            if r == radius 
               arrayCirc(1,count,:) = temp(y,x,:);
               count = count + 1;
            end
        end
    end
        
    % Values in a circle that overlaps an OBSTACLE or a region OUTSIDE the monolayer are set to 0.  
    [rowNew, columnNew, frameNew] = size(arrayCirc); 
    random = zeros(rowNew,columnNew,frameNew); 
    for a = 1:columnNew
        if (arrayCirc(1,a,t) == arrayCirc(1,a,t+1)) && (arrayCirc(1,a,t+1) == arrayCirc(1,a,t+2))
           arrayCirc = random;
        end 
    end
    
    %% Calculate the phase difference, sum of phase difference, and winding number.
    sumPhase = 0; 
    [~,n,~] = size(arrayCirc); 
    difPhase = zeros(1,n);

    for i = 1:n-1 
        difPhase(1,i) = arrayCirc(1, i+1, t) - arrayCirc(1, i, t); 

        % A large difPhase value is considered to be due to a crossing
        % over of the 0 <-> 360 axis. So 360 is added or subtracted. -> I
        % think this section might be bringing about some false positives. 
        % For example, in recordings that are not so clear, it may be that
        % difPhase really is greater than 180 or -180, but is being
        % inaccurately compensated. 
        if difPhase(1,i) < -180 
           difPhase(1,i) = difPhase(1,i) + 360;
        elseif difPhase(1,i) > 180
           difPhase(1,i) = difPhase(1,i) - 360;
        end 

        sumPhase = sumPhase + difPhase(1,i);
    end

    % Diff btw the 1st and last pixels 
    difPhase(1,n) = arrayCirc(1, 1, t) - arrayCirc(1, n, t); 

    if difPhase(1,n) < -180
       difPhase(1,n) = difPhase(1,n) + 360;
       elseif difPhase(1,n) > 180
       difPhase(1,n) = difPhase(1,n) - 360;
    end 

    sumPhase = sumPhase + difPhase(1,n); 
    windingNum = sumPhase/360;

    % Save the values for which winding number = -1 or +1 so it can be mapped.
        if windingNum < -0.99 && windingNum > - 1.01            
           clock(1,clocktic) = centerX;
           clock(2,clocktic) = centerY;
           clocktic = clocktic+1;                
           elseif windingNum < 1.01 && windingNum > 0.99  
           counterClock(1,countertic) = centerX; 
           counterClock(2,countertic) = centerY; 
           countertic = countertic+1;                
        end        
    end
end 
 
%% Map the points that give an windingNum of -1 or +1 
    figure(3);
    subplot(1,12,2:12);  % This subplots the monolayer on the right side of the bar that indicates the frames 
    pcolor(video(:,:,t));  
    colormap(prism(5)) % Max number of colors for prism is 6!     
    shading flat 
    hold on

    set(gca, 'YTick', 0:10:80); 
    set(gca,'YDir','reverse'); % This reverses the direction of counting on Y-axis 
    set(gca, 'XTick', 0:10:80); 
    axis square; % This makes the monolayer into a square 
    grid on
    grid minor
    set(gca,'layer','top');
    xlabel('Pixel'); ylabel('Pixel');
    
    for i = 1:speed:size(clock,2) % If there's nothing in clock, it returns an error saying 'data exceeds matrix dimensions' 
        if clock(1,i) == 0
            continue 
        end
        x_px = clock(1,i); 
        y_px = clock(2,i);
        plot(x_px,y_px,'k*');
        sourceMap(y_px,x_px, 1) = sourceMap(y_px,x_px, 1) - 1;
        hold on
    end

    for i = 1:speed:size(counterClock,2) % If there's nothing in counterClock, it returns an error saying 'data exceeds matrix dimensions' 
        if counterClock(1,i) == 0
            continue 
        end
        x_px = counterClock(1,i); 
        y_px = counterClock(2,i);
        plot(x_px,y_px,'w*');
        sourceMap(y_px,x_px, 2) = sourceMap(y_px,x_px, 2) + 1;
        hold on 
    end

    % Counter for cycle
    if t == finFrame %iniFrame + (cycle*ceil(Period*Fs))
        x_px = column/10;
        y_px = row/10;
        plot(x_px,y_px,'o');
        hold on
        cycle = cycle + 1;
    end
    hold off
    
    % Bar 
    subplot(1,12,1);
    plot(1,t,'g.','MarkerSize',20);
    set(gca,'XTick',NaN);
    ylim([iniFrame finFrame]); % Sets y-axis tick marks 
    ylabel('Frame Number');
    hold on;
    
    % The function that creates the video 
    F(t)=getframe(gcf);

    % Saves .gif file 
    if save==1;
        drawnow
        im = frame2im(F(t));
        [A,map] = rgb2ind(im,256); 
        
        if t == iniFrame;
            imwrite(A,map,'testing_gif.gif','gif','LoopCount',Inf,'DelayTime',T*speed*10); % multiplied by 10 to slow it down. 
            else
            imwrite(A,map,'testing_gif.gif','gif','WriteMode','append','DelayTime',T*speed*10);
        end
    end
end

%% Source map 
min = -1; 
max = 1; 
NewSourceMap = sourceMap/((finFrame-(iniFrame-1))/speed); 
figure(1);
imagesc(NewSourceMap(:,:,1)); 
axis square; 
colormap(hot);
caxis([min, max])
colorbar;

figure(2); 
imagesc(NewSourceMap(:,:,2)); 
axis square;
colormap(hot);
caxis([min, max])
colorbar;

toc