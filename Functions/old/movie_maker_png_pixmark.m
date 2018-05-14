function [pic]=movie_maker_png_pixmark(g,framenum1,framenum2,ff, pix1, pix2, pix3, pix4)

% framenum1,2 are the frames you want to start and stop the movie at.
% ff is the number of frames you skip
%pix1,2,3,4 are pixels of interest that are marked on the video.

close all;

g(pix1(2), pix1(1),:)=NaN; 
g(pix2(2), pix2(1),:)=NaN; 
g(pix3(2), pix3(1),:)=NaN; 
g(pix4(2), pix4(1),:)=NaN; 

frame=framenum1:ff:framenum2;

for i=1:length(frame)
    hold on
    subplot(1,12,1);   
    plot(1,frame(i),'g.','MarkerSize',20);
    xlim([0 2]); ylim([framenum1 framenum2]);
    ylabel('Frame Number');
    subplot(1,12,2:12);
    pcolor(g(:,:,frame(i)));
    shading flat
    set(gca,'YDir','reverse');
    axis square;colormap(summer) 
    xlabel('Pixel'); ylabel('Pixel');
    cmin=min(g(20,40,:));
    cmax=max(g(60,40,:));
    caxis([cmin cmax]);
    pic=getframe(gcf);
    imwrite(pic.cdata,sprintf('~/Videos/%04d.png',i));
    hold off;
end

%time1=framenum1/40-1/40;
%time2=framenum2/40-1/40;
