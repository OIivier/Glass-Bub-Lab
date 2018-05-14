function movie_maker_mod(g,framenum1,framenum2,speed,x_px,y_px,x_px2,y_px2,IC,switchtime,save)

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

% if pix1(1)~=0;
% g(pix1(2), pix1(1),:)=NaN; 
% g(pix2(2), pix2(1),:)=NaN; 
% g(pix3(2), pix3(1),:)=NaN; 
% g(pix4(2), pix4(1),:)=NaN;
% end


% Section 1: Color the map 
for i=1:length(frame)
    figure(6);

    imagesc(g(:,:,frame(i)));  
    %pcolor colors the matrix based on default colors. The default shading 
    %is faceted, and this will delete the last row and column of C.  
    %imagesc colors the m-by-n grid of pixels according to the colormap. 
    %This function will keep all of the m-by-n pixels. 
    %shading flat;
    colormap(summer) % Another option is Hot
    cmin=min(g(35,50,:)); % Find the minimum at this point over time
    cmax=max(g(35,50,:)); % Find the maximum at this point over time 
%     cmin=max(min(g(20,30,:)),mean(g(20,30,:))-1*std(g(20,30,:)));
%     cmax=min(max(g(20,30,:)),mean(g(20,30,:))+1*std(g(20,30,:)));
    caxis([cmin cmax]);  % Set the min and max colors to be based on found min & max
%    caxis([-1.97 1.62]);

    %Notes on Shading 
    %flat: each mesh line segment and face has a constant color determined 
    %by the color value at the endpoint of the segment or the corner of the 
    %face that has the smallest index or indices.
    %interp: resamples data at a higher rate using lowpass interpolation. 
    %Y = interp(X,R) resamples the sequence in vector X at R times the 
    %original sample rate. With shading interp, each cell is colored by 
    %bilinear interpolation of the colors at its four vertices, using all 
    %elements of C.
    
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
    
    if x_px ~= 0
        plot(x_px,y_px,'wo');
    end
    
    if x_px2 ~= 0
        plot(x_px2,y_px2,'wo');
    end

%    plot([35 55], [68 55], '-o')                       % <--
%    plot(55,55, 'wo')                       % <--
%    plot(55,40, 'w*')                       % <--

%    text(2.5,1,[sprintf('%.3f',frame(i)*0.025) 's'],      ...   % <--
    

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
           imwrite(A,map,'/home/vincent/Documents/Research Thesis/FiguresAndGifs/testing_gif.gif','gif','LoopCount',Inf,'DelayTime',T*speed);
           else
           imwrite(A,map,'/home/vincent/Documents/Research Thesis/FiguresAndGifs/testing_gif.gif','gif','WriteMode','append','DelayTime',T*speed);
        end
    end
        
    
%    if save==1;
%       movie2avi(F,'inputE2','fps',40);    
%    end
end

toc
