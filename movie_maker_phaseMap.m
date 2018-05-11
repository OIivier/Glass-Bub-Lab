function movie_maker_phaseMap(g,framenum1,framenum2,speed,x_px,y_px,x_px1,y_px1,IC,switchtime,save)

% framenum1,2 are the frames you want to start and stop the movie at
% speed=2, and fps=20 let you make an avi file (filename=input.avi) that is
% in realtime

% IC=1 means that initially the condition is spontaneous with no perfusion
% IC=2 means that initially the condition is spontaneous with perfusion

% If switchtime~=0 then you will get a change in color after the time (in seconds)

% If save==1 the movie will be saved, otherwise you can just view the movie

Fs = 40;              % Sampling frequency
T = 1/Fs;             % Sampling period

tic
switchtime=switchtime*Fs;
close all;

frame = framenum1:speed:framenum2; %speed determines how many frames are skipped

% Section 1: Color the map 

for i=1:length(frame)
    figure(6);
    pcolor(g(:,:,frame(i)));
%    imagesc(g(:,:,frame(i)));
    C=colormap(prism(5)); % Max number of colors for prism is 6! 
    caxis([-180, 180])
    shading flat
    
%     caxis([-180, 324])
%     C = [C; 0 0 0; 1 1 1]; % 1 1 1: white, 0 0 0: black
    colormap(C);
    %For a custom colormap, specify map as a three-column matrix of values 
    %in the range [0,1] where each row is an RGB triplet that defines one 
    %color of the colormap. An RGB triplet is a three-element row vector 
    %specifying the red, green, and blue intensities for a color.
    
% Section 2: Draws the bar  
    subplot(1,12,1);
    
    if switchtime~=0
        if IC==1
            if frame(i)<=switchtime
                plot(1,frame(i),'g.','MarkerSize',20);
            end
            if frame(i)>switchtime
                 plot(1,frame(i),'b.','MarkerSize',20);
            end
        end
        if IC==2
            if frame(i)<=switchtime
                plot(1,frame(i),'b.','MarkerSize',20);
            end
            if frame(i)>switchtime
                 plot(1,frame(i),'r.','MarkerSize',20);
            end
        end
    end
    
    if switchtime==0;
        if IC==1
        plot(1,frame(i),'g.','MarkerSize',20);
        end
        if IC==2
        plot(1,frame(i),'b.','MarkerSize',20);
        end
        if IC==3
        plot(1,frame(i),'r.','MarkerSize',20);
        end
    end

    set(gca,'XTick',NaN);
    ylim([framenum1 framenum2]); % Sets y-axis tick marks 
    ylabel('Frame Number');
    hold on;
    
% Section 3: Plots the Monolayer to the right of the Frame bar 
    subplot(1,12,2:12);  % This subplots the monolayer on the right side of the bar that indicates the frames 
    hold on;
    set(gca, 'YTick', 0:10:80); 
    set(gca,'YDir','reverse'); % This reverses the direction of counting on Y-axis 
    set(gca, 'XTick', 0:10:80); 
    axis square; % This makes the box that the monolayer is in a square 
    grid on
    grid minor
    set(gca,'layer','top');
    xlabel('Pixel'); ylabel('Pixel');
    
%     if x_px ~= 0
%         plot(x_px,y_px,'w*');
%     end

    if x_px ~= 0
        plot(x_px,y_px,'ko', x_px1,y_px1,'ko', x_px2,y_px2,'ko', x_px3,y_px3,'ko');
    end

%     text(2.5,1,[sprintf('%.3f',frame(i)*0.025) 's'],      ...   % <--
      text(77,2,[sprintf('%.3f',frame(i)*0.025) 's'],      ...   % <--
        'Color', 'w',                 ...
        'FontSize', 12,               ...
        'HorizontalAlignment','right',...   
        'VerticalAlignment','cap');         
    
    hold off;
    F(i)=getframe(gcf);
%   waitforbuttonpress                      % <--

% Section 4: Creates .gif file 
    if save==1;
        drawnow
        im = frame2im(F(i));
        [A,map] = rgb2ind(im,256); 
       
        if i == 1;
            imwrite(A,map,'testing_gif.gif','gif','LoopCount',Inf,'DelayTime',T*speed);
        else
            imwrite(A,map,'testing_gif.gif','gif','WriteMode','append','DelayTime',T*speed);
        end
    end
end

toc
