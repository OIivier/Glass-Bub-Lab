function [pic]=movie_maker_png(g,framenum1,framenum2,ff, pix1, pix2, pix3, pix4)

% framenum1,2 are the frames you want to start and stop the movie at

close all;

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
    cmax=max(g(20,40,:));
    caxis([cmin cmax]);
    pic=getframe(gcf);
    imwrite(pic.cdata,sprintf('~/Videos/%04d.png',i));
    hold off;
end

%time1=framenum1/40-1/40;
%time2=framenum2/40-1/40;
