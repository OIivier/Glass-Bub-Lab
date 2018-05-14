%IBI analysis script started April 2009

%should close all opened figures before running this, but jsut in case make
%a fresh one
figure

%load file
%[B,A,M,N,NumFrames,NumCols,NumRows,FrameInt,AcquisitionRatio,ibis,ActMap,Times] = daRead4(daFile,XCoord,YCoord);
%from plot determine where to set threshold for IBIs and put it in here
thresh=7;
threshdoc=0.45;
%truncate first n frames
n=1;
N2=SF1(:,:,n:end);

c
%calculate and plot timeseries and IBIs for some pixels
[IBI1,IBIdoc1,time1,timedoc1,meanIBI1,meanIBIdoc1,stdIBI1,stdIBIdoc1]=ibiBBplot1(N2,27,31,thresh,threshdoc);
[IBI2,IBIdoc2,time2,timedoc2,meanIBI2,meanIBIdoc2,stdIBI2,stdIBIdoc2]=ibiBBplot1(N2,32,51,thresh,threshdoc);
[IBI3,IBIdoc3,time3,timedoc3,meanIBI3,meanIBIdoc3,stdIBI3,stdIBIdoc3]=ibiBBplot1(N2,55,59,thresh,threshdoc);
[IBI4,IBIdoc4,time4,timedoc4,meanIBI4,meanIBIdoc4,stdIBI4,stdIBIdoc4]=ibiBBplot1(N2,64,32,thresh,threshdoc);
%[IBI5,IBIdoc5,time5,timedoc5,meanIBI5,meanIBIdoc5,stdIBI5,stdIBI5]=ibiBBplot1(N2,43,31,thresh,threshdoc);
%[IBI6,IBIdoc6,time6,timedoc6,meanIBI6,meanIBIdoc6,stdIBI6,stdIBIdoc6]=ibiBBplot1(N2,44,31,thresh,threshdoc);
%[IBI7,IBIdoc7,time7,timedoc7,meanIBI7,meanIBIdoc7,stdIBI7,stdIBIdoc7]=ibiBBplot1(N2,47,38,thresh,threshdoc);
%[IBI8,IBIdoc8,time8,timedoc8,meanIBI8,meanIBIdoc8,stdIBI8,stdIBIdoc8]=ibiBBplot1(N2,48,40,thresh,threshdoc);
%[IBI9,IBIdoc9,time9,timedoc9,meanIBI9,meanIBIdoc9,stdIBI9,stdIBIdoc9]=ibiBBplot1(N1,60,60,thresh,threshdoc);

%IBImeanvec=zeros(1,1681);
%IBIstdvec=zeros(1,1681);
%IBImeanspace=zeros(41,41);
%IBImeanspace=zeros(41,41);

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