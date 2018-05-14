%this function plots all IBIs from the different conditions
%for a particular case (the first one I built) see allIBIplot.m
%****remember to name your q so as to remember pixel. ex.: allIBI2060
%(i couldn't make the varname a workspace variable) xD)

%% pix1 %%


[SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/spont.da');
[IBI1, IBIdoc1, time1, timedoc1, meanIBI1, meanIBIdoc1, stdIBI1, stdIBIdoc1, IBI2, IBIdoc2, time2, timedoc2, meanIBI2, meanIBIdoc2, stdIBI2, stdIBIdoc2, IBI3, IBIdoc3, time3, timedoc3, meanIBI3, meanIBIdoc3, stdIBI3, stdIBIdoc3, IBI4, IBIdoc4, time4, timedoc4, meanIBI4, meanIBIdoc4, stdIBI4, stdIBIdoc4]=IBI_multipixel_analysis2(SF1, [30 60], [20 20], [20 65], [20 18], 30, 30, 30, 30, 0.4);

allIBIpix1s=[IBIdoc1;];
allIBIpix2s=[IBIdoc2;];
allIBIpix3s=[IBIdoc3;];
allIBIpix4s=[IBIdoc4;];


for c=2:3;
    if c==2;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/spont01.da');
    end
    if c==3;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/spont02.da');
    end
    [IBI1, IBIdoc1, time1, timedoc1, meanIBI1, meanIBIdoc1, stdIBI1, stdIBIdoc1, IBI2, IBIdoc2, time2, timedoc2, meanIBI2, meanIBIdoc2, stdIBI2, stdIBIdoc2, IBI3, IBIdoc3, time3, timedoc3, meanIBI3, meanIBIdoc3, stdIBI3, stdIBIdoc3, IBI4, IBIdoc4, time4, timedoc4, meanIBI4, meanIBIdoc4, stdIBI4, stdIBIdoc4]=IBI_multipixel_analysis2(SF1, [30 60], [20 20], [30 35], [20 18], 30, 30, 30, 30, 0.4);
    allIBIpix1s=[allIBIpix1s; IBIdoc1;];
    allIBIpix2s=[allIBIpix2s; IBIdoc2;];
    allIBIpix3s=[allIBIpix3s; IBIdoc3;];
    allIBIpix4s=[allIBIpix4s; IBIdoc4;];
end

allIBIpix1c=[];
allIBIpix2c=[];
allIBIpix3c=[];
allIBIpix4c=[];

for c=4:6
    if c==4;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/control.da');
    end
    if c==5;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/control01.da');
    end
    if c==6;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/control02.da');
    end
    [IBI1, IBIdoc1, time1, timedoc1, meanIBI1, meanIBIdoc1, stdIBI1, stdIBIdoc1, IBI2, IBIdoc2, time2, timedoc2, meanIBI2, meanIBIdoc2, stdIBI2, stdIBIdoc2, IBI3, IBIdoc3, time3, timedoc3, meanIBI3, meanIBIdoc3, stdIBI3, stdIBIdoc3, IBI4, IBIdoc4, time4, timedoc4, meanIBI4, meanIBIdoc4, stdIBI4, stdIBIdoc4]=IBI_multipixel_analysis2(SF1, [30 60], [20 20], [30 35], [20 18], 30, 30, 30, 30, 0.4);
    allIBIpix1c=[allIBIpix1c; IBIdoc1;];
    allIBIpix2c=[allIBIpix2c; IBIdoc2;];
    allIBIpix3c=[allIBIpix3c; IBIdoc3;];
    allIBIpix4c=[allIBIpix4c; IBIdoc4;];
end


allIBIpix1d=[];
allIBIpix2d=[];
allIBIpix3d=[];
allIBIpix4d=[];

for c=7:25
    if c==7;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol.da');
        threshs = [25 25 25 25];
    end
    if c==8;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol01.da');
        threshs = [30 30 40 30];
    end
    if c==9;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol02.da');
        threshs = [50 50 50 50];
    end
    if c==10;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/Dish4/E40310p75umol03.da');
        threshs = [30 30 30 30];
    end
    if c==11;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol04.da');
        threshs = [30 30 30 30];
    end
    if c==12;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol06.da');
        threshs = [30 30 30 30];
    end
    if c==13;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol07.da');
        threshs = [30 30 30 30];
    end
    if c==14;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol08.da');
        threshs = [30 30 30 30];
    end
    if c==15;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol09.da');
        threshs = [30 30 30 30];
    end
    if c==16;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol10.da');
        threshs = [30 30 30 30];
    end
    if c==17;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol11.da');
        threshs = [30 30 30 30];
    end
    if c==18;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol12.da');
        threshs = [30 30 30 30];
    end
    if c==19;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol13.da');
        threshs = [30 30 30 30];
    end
    if c==20;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/Dish4/E40310p75umol14.da');
        threshs = [30 30 30 30];
    end
    if c==21;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol15.da');
        threshs = [30 30 30 30];
    end
    if c==22;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol16.da');
        threshs = [30 30 30 30];
    end
    if c==23;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol17.da');
        threshs = [30 30 30 30];
    end
    if c==24;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol18.da');
        threshs = [30 30 30 30];
    end
    if c==25;
        [SF1] = daRead7('/cnd0/home/borek/expt/june24/Dish3/E40310p75umol19.da');
        threshs = [30 30 30 30];
    end
    if c==26;
        [SF1] = daRead7('/cnd0/home/gabriels/Desktop/Experiments/june17/Dish4/E40310p75umol20.da');
        threshs = [30 30 30 30];
    end
    if c==27;
        [SF1] = daRead7('/cnd0/home/gabriels/Desktop/Experiments/june17/Dish4/E40310p75umol21.da');
        threshs = [65 65 65 65];
    end
    
    
    [IBI1, IBIdoc1, time1, timedoc1, meanIBI1, meanIBIdoc1, stdIBI1, stdIBIdoc1, IBI2, IBIdoc2, time2, timedoc2, meanIBI2, meanIBIdoc2, stdIBI2, stdIBIdoc2, IBI3, IBIdoc3, time3, timedoc3, meanIBI3, meanIBIdoc3, stdIBI3, stdIBIdoc3, IBI4, IBIdoc4, time4, timedoc4, meanIBI4, meanIBIdoc4, stdIBI4, stdIBIdoc4]=IBI_multipixel_analysis2(SF1, [30 60], [20 20], [30 35], [20 18], threshs(1), threshs(2), threshs(3), threshs(4), 0.4);
   
    allIBIpix1d=[allIBIpix1d; IBIdoc1;];
    allIBIpix2d=[allIBIpix2d; IBIdoc2;];
    allIBIpix3d=[allIBIpix3d; IBIdoc3;];
    allIBIpix4d=[allIBIpix4d; IBIdoc4;];
end

allIBIpix1w=[];
allIBIpix2w=[];
allIBIpix3w=[];
allIBIpix4w=[];




spont1=1:length(allIBIpix1s);
spont2=1:length(allIBIpix2s);
spont3=1:length(allIBIpix3s);
spont4=1:length(allIBIpix4s);

control1=1:length(allIBIpix1c);
control2=1:length(allIBIpix2c);
control3=1:length(allIBIpix3c);
control4=1:length(allIBIpix4c);

drug1=1:length(allIBIpix1d);
drug2=1:length(allIBIpix2d);
drug3=1:length(allIBIpix3d);
drug4=1:length(allIBIpix4d);

wash1=1:length(allIBIpix1w);
wash2=1:length(allIBIpix2w);
wash3=1:length(allIBIpix3w);
wash4=1:length(allIBIpix4w);

figure(1);
subplot(4,1,1);
plot(spont1,allIBIpix1s,'g*');
hold on;
plot((length(spont1)+1:(length(spont1)+length(control1))),allIBIpix1c,'b*');
plot(((length(spont1)+length(control1)+1):(length(spont1)+length(control1)+length(drug1))),allIBIpix1d,'r*');
plot(((length(spont1)+length(control1)+length(drug1)+1):(length(spont1)+length(control1)+length(drug1)+length(wash1))),allIBIpix1w,'c*');
title('All IBIs June 17 Dish 4 Pix 1','FontSize',30); ylabel('IBI(s)','FontSize',16);

subplot(4,1,2);
plot(spont2,allIBIpix2s,'g*');
hold on;
plot((length(spont2)+1:(length(spont2)+length(control2))),allIBIpix2c,'b*');
plot(((length(spont2)+length(control2)+1):(length(spont2)+length(control2)+length(drug2))),allIBIpix2d,'r*');
plot(((length(spont2)+length(control2)+length(drug2)+1):(length(spont2)+length(control2)+length(drug2)+length(wash2))),allIBIpix2w,'c*');
title('All IBIs June 17 Dish 4 Pix 2','FontSize',30); ylabel('IBI(s)','FontSize',16);

subplot(4,1,3);
plot(spont3,allIBIpix3s,'g*');
hold on;
plot((length(spont3)+1:(length(spont3)+length(control3))),allIBIpix3c,'b*');
plot(((length(spont3)+length(control3)+1):(length(spont3)+length(control3)+length(drug3))),allIBIpix3d,'r*');
plot(((length(spont3)+length(control3)+length(drug3)+1):(length(spont3)+length(control3)+length(drug3)+length(wash3))),allIBIpix3w,'c*');
title('All IBIs June 17 Dish 4 Pix 3','FontSize',30); ylabel('IBI(s)','FontSize',16);

subplot(4,1,4);
plot(spont4,allIBIpix4s,'g*');
hold on;
plot((length(spont4)+1:(length(spont4)+length(control4))),allIBIpix4c,'b*');
plot(((length(spont4)+length(control4)+1):(length(spont4)+length(control4)+length(drug4))),allIBIpix4d,'r*');
plot(((length(spont4)+length(control4)+length(drug4)+1):(length(spont4)+length(control4)+length(drug4)+length(wash4))),allIBIpix4w,'c*');
title('All IBIs June 17 Dish 4 Pix 4','FontSize',30); ylabel('IBI(s)','FontSize',16); xlabel('Beat Number Analyzed','FontSize',16);


%% barograms %%

meanp1s=mean(allIBIpix1s); stanp1s=std(allIBIpix1s);
meanp2s=mean(allIBIpix2s); stanp2s=std(allIBIpix2s);
meanp3s=mean(allIBIpix3s); stanp3s=std(allIBIpix3s);
meanp4s=mean(allIBIpix4s); stanp4s=std(allIBIpix4s);

meanp1c=mean(allIBIpix1c); stanp1c=std(allIBIpix1c);
meanp2c=mean(allIBIpix2c); stanp2c=std(allIBIpix2c); 
meanp3c=mean(allIBIpix3c); stanp3c=std(allIBIpix3c);
meanp4c=mean(allIBIpix4c); stanp4c=std(allIBIpix4c);

meanp1d=mean(allIBIpix1d); stanp1d=std(allIBIpix1d);
meanp2d=mean(allIBIpix2d); stanp2d=std(allIBIpix2d); 
meanp3d=mean(allIBIpix3d); stanp3d=std(allIBIpix3d);
meanp4d=mean(allIBIpix4d); stanp4d=std(allIBIpix4d);

meanp1w=mean(allIBIpix1w); stanp1w=std(allIBIpix1w);
meanp2w=mean(allIBIpix2w); stanp2w=std(allIBIpix2w);
meanp3w=mean(allIBIpix3w); stanp3w=std(allIBIpix3w);
meanp4w=mean(allIBIpix4w); stanp4w=std(allIBIpix4w);

numbins=10;
figure(2);
subplot(2,2,1);
bar(allIBIpix1s,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','g','EdgeColor','w');
xlim([0 4]);
title(['Time Series for Pixel 1 (mean,std) = (',num2str(meanp1s),',',num2str(stanp1s),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,2);
bar(allIBIpix1c,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');
xlim([0 4]);
title(['Time Series for Pixel 1 (mean,std) = (',num2str(meanp1c),',',num2str(stanp1c),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,3);
bar(allIBIpix1d,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
xlim([0 4]);
title(['Time Series for Pixel 1 (mean,std) = (',num2str(meanp1d),',',num2str(stanp1d),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,4);
bar(allIBIpix1w,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w');
xlim([0 4]);
title(['Time Series for Pixel 1 (mean,std) = (',num2str(meanp1w),',',num2str(stanp1w),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);

figure(3);
subplot(2,2,1);
bar(allIBIpix2s,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','g','EdgeColor','w');
title(['Time Series for Pixel 2 (mean,std) = (',num2str(meanp2s),',',num2str(stanp2s),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,2);
bar(allIBIpix2c,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');
title(['Time Series for Pixel 2 (mean,std) = (',num2str(meanp2c),',',num2str(stanp2c),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,3);
bar(allIBIpix2d,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
title(['Time Series for Pixel 2 (mean,std) = (',num2str(meanp2d),',',num2str(stanp2d),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,4);
bar(allIBIpix2w,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w');
title(['Time Series for Pixel 2 (mean,std) = (',num2str(meanp2w),',',num2str(stanp2w),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);

figure(4);
subplot(2,2,1);
bar(allIBIpix3s,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','g','EdgeColor','w');
title(['Time Series for Pixel 3 (mean,std) = (',num2str(meanp3s),',',num2str(stanp3s),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,2);
bar(allIBIpix3c,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');
title(['Time Series for Pixel 3 (mean,std) = (',num2str(meanp3c),',',num2str(stanp3c),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,3);
bar(allIBIpix3d,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
title(['Time Series for Pixel 3 (mean,std) = (',num2str(meanp3d),',',num2str(stanp3d),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,4);
bar(allIBIpix3w,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w');
title(['Time Series for Pixel 3 (mean,std) = (',num2str(meanp3w),',',num2str(stanp3w),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);

figure(5);
subplot(2,2,1);
bar(allIBIpix4s,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','g','EdgeColor','w');
title(['Time Series for Pixel 4 (mean,std) = (',num2str(meanp4s),',',num2str(stanp4s),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,2);
bar(allIBIpix4c,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');
title(['Time Series for Pixel 4 (mean,std) = (',num2str(meanp4c),',',num2str(stanp4c),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,3);
bar(allIBIpix4d,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
title(['Time Series for Pixel 4 (mean,std) = (',num2str(meanp4d),',',num2str(stanp4d),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
subplot(2,2,4);
bar(allIBIpix4w,numbins);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','c','EdgeColor','w');
title(['Time Series for Pixel 4 (mean,std) = (',num2str(meanp4w),',',num2str(stanp4w),')'],'FontSize',20,'Interpreter','none');
ylabel('Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);


[n1,xout1]=bar(allIBIpix4s,numbins); n1=n1./length(allIBIpix4s);
[n2,xout2]=bar(allIBIpix4c,numbins); n2=n2./length(allIBIpix4c);
[n3,xout3]=bar(allIBIpix4d,numbins); n3=n3./length(allIBIpix4d);
[n4,xout4]=bar(allIBIpix4w,200); n4=n4./length(allIBIpix4w);
figure(7);
stem(xout1,n1,'g');
hold on;
stem(xout2,n2,'b');
stem(xout3,n3,'r');
stem(xout4,n4,'c');
ylabel('Normalized Number of Observations','FontSize',12); xlabel('IBI (s)','FontSize',12);
legend('Spontaneous','Control','E4031 (0.75 umol)','Wash');