function movie_maker(g,framenum1,framenum2,speed,pix1,pix2,pix3,pix4,IC,switchtime,save)

% framenum1,2 are the frames you want to start and stop the movie at
% speed=2, and fps=20 let you make an avi file (filename=input.avi) that is
% in realtime

% pix1-4 are pixels of interest to be whited-out
% If pix1(1)=0 then none of the pixels will be whited out (you still need to put values in for the other pixels though)

% IC=1 means that initially the condition is spontaneous with no perfusion
% IC=2 means that initially the condition is spontaneous with perfusion


% If switchtime~=0 then you will get a change in color after the time (in seconds)

% If save==1 the movie will be saved, otherwise you can just view the movie

tic
switchtime=switchtime*40;
close all;

frame=framenum1:speed:framenum2;

if pix1(1)~=0;
g(pix1(2), pix1(1),:)=NaN; 
g(pix2(2), pix2(1),:)=NaN; 
g(pix3(2), pix3(1),:)=NaN; 
g(pix4(2), pix4(1),:)=NaN;
end

for i=1:length(frame)
    figure(6);
    pcolor(g(:,:,frame(i)));
    shading interp
    colormap(summer) % Another option is Hot
    % cmin=min(min(min(g))); %these make limits too big
    % cmax=max(max(max(g)));
    cmin=min(g(20,41,:));
    cmax=max(g(20,41,:));    
    caxis([cmin cmax])
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
    
    set(gca,'XTick',NaN); ylim([framenum1 framenum2]);
    ylabel('Frame Number');
    hold on;
    subplot(1,12,2:12);
    hold on;
    set(gca,'YDir','reverse');
    axis square;
    xlabel('Pixel'); ylabel('Pixel');
    hold off;
    F(i)=getframe(gcf);
end

if save==1;
movie2avi(F,'inputE2','fps',40);
end

toc
