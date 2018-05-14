%IBI analysis script started April 2009

%should close all opened figures before running this, but jsut in case make
%a fresh one
figure

%load file
%[B,A,M,N,NumFrames,NumCols,NumRows,FrameInt,AcquisitionRatio,ibis,ActMap,Times] = daRead4(daFile,XCoord,YCoord);
%from plot determine where to set threshold for IBIs and put it in here
thresh=10;
threshdoc=0.55;
%truncate first 30seconds
N1=N(:,:,600:end);




%calculate and plot timeseries and IBIs for some pixels
[IBI2020,IBIdoc2020,time2020,timedoc2020,meanIBI2020,meanIBIdoc2020,stdIBI2020,stdIBIdoc2020]=ibiBBplot(N1,20,20,thresh,threshdoc);
[IBI2040,IBIdoc2040,time2040,timedoc2040,meanIBI2040,meanIBIdoc2040,stdIBI2040,stdIBIdoc2040]=ibiBBplot(N1,20,40,thresh,threshdoc);
[IBI2060,IBIdoc2060,time2060,timedoc2060,meanIBI2060,meanIBIdoc2060,stdIBI2060,stdIBIdoc2060]=ibiBBplot(N1,20,60,thresh,threshdoc);
[IBI4020,IBIdoc4020,time4020,timedoc4020,meanIBI4020,meanIBIdoc4020,stdIBI4020,stdIBIdoc4020]=ibiBBplot(N1,40,20,thresh,threshdoc);
[IBI4040,IBIdoc4040,time4040,timedoc4040,meanIBI4040,meanIBIdoc4040,stdIBI4040,stdIBIdoc4040]=ibiBBplot(N1,40,40,thresh,threshdoc);
[IBI4060,IBIdoc4060,time4060,timedoc4060,meanIBI4060,meanIBIdoc4060,stdIBI4060,stdIBIdoc4060]=ibiBBplot(N1,40,60,thresh,threshdoc);
[IBI6020,IBIdoc6020,time6020,timedoc6020,meanIBI6020,meanIBIdoc6020,stdIBI6020,stdIBIdoc6020]=ibiBBplot(N1,60,20,thresh,threshdoc);
[IBI6040,IBIdoc6040,time6040,timedoc6040,meanIBI6040,meanIBIdoc6040,stdIBI6040,stdIBIdoc6040]=ibiBBplot(N1,60,40,thresh,threshdoc);
[IBI6060,IBIdoc6060,time6060,timedoc6060,meanIBI6060,meanIBIdoc6060,stdIBI6060,stdIBIdoc6060]=ibiBBplot(N1,60,60,thresh,threshdoc);

IBImeanvec=zeros(1,1681);
IBIstdvec=zeros(1,1681);
IBImeanspace=zeros(41,41);
IBImeanspace=zeros(41,41);

%for y=20:60;
%    for x=20:60;
%        [IBIdoc,timedoc,meanIBIdoc,stdIBIdoc]=ibiBB(N1,x,y,thresh,threshdoc);
%        if meanIBIdoc<
%        IBImeanspace(x-19,y-19)=meanIBIdoc;
%        IBIstdspace(x-19,y-19)=stdIBIdoc;
%        IBImeanvec((x-20+1)+41*(y-20))=meanIBIdoc;
%        IBIstdvec((x-20+1)+41*(y-20))=stdIBIdoc;
%    end
%end

%figure
%hist(IBImeanspace)
%hist(IBIstdspace)
        