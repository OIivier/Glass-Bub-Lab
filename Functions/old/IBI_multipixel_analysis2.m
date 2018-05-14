function [IBI1, IBIdoc1, time1, timedoc1, meanIBI1, meanIBIdoc1, stdIBI1, stdIBIdoc1, IBI2, IBIdoc2, time2, timedoc2, meanIBI2, meanIBIdoc2, stdIBI2, stdIBIdoc2, IBI3, IBIdoc3, time3, timedoc3, meanIBI3, meanIBIdoc3, stdIBI3, stdIBIdoc3, IBI4, IBIdoc4, time4, timedoc4, meanIBI4, meanIBIdoc4, stdIBI4, stdIBIdoc4]=IBI_multipixel_analysis2(SF1, pix1, pix2, pix3, pix4, threshspix1, threshspix2, threshspix3, threshspix4, threshdoc)

close all;

%IBI analysis script started April 2009

%load file
%[B,A,M,N,NumFrames,NumCols,NumRows,FrameInt,AcquisitionRatio,ibis,ActMap,Times] = daRead4(daFile,XCoord,YCoord);
%from plot determine where to set threshold for IBIs and put it in here

%calculate and plot timeseries and IBIs for some pixels
[IBI1,IBIdoc1,time1,timedoc1,meanIBI1,meanIBIdoc1,stdIBI1,stdIBIdoc1]=ibiBBplot2(SF1,pix1(1),pix1(2),threshspix1,threshdoc);
[IBI2,IBIdoc2,time2,timedoc2,meanIBI2,meanIBIdoc2,stdIBI2,stdIBIdoc2]=ibiBBplot2(SF1,pix2(1),pix2(2),threshspix2,threshdoc);
[IBI3,IBIdoc3,time3,timedoc3,meanIBI3,meanIBIdoc3,stdIBI3,stdIBIdoc3]=ibiBBplot2(SF1,pix3(1),pix3(2),threshspix3,threshdoc);
[IBI4,IBIdoc4,time4,timedoc4,meanIBI4,meanIBIdoc4,stdIBI4,stdIBIdoc4]=ibiBBplot2(SF1,pix4(1),pix4(2),threshspix4,threshdoc);

end