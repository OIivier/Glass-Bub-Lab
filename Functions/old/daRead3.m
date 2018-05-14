function [B,A,C,M,M1,N,N1,NumFrames,NumCols,NumRows,FrameInt,AcquisitionRatio,ActMap] = daRead3(daFile,XCoord,YCoord)
%this function plots the time series for a single pixel of the redshirt
%camera, both filtered and unfiltered. Pass it the file name and 
%coordinates. (1,1) is top left. It returns three matrices:M, N, ibis
%representing unfiltered data, filtered data, and a binary 3D matrix
%of size M where a 1 denotes an upstroke
tic
data=fopen(daFile);
file=fread(data,'bit16');
NumFrames=file(5)
NumCols=file(385)
NumRows=file(386)
FrameInt=file(389)/1000000
AcquisitionRatio=file(392)
INDEX = 1:1:NumFrames;
Times = INDEX*FrameInt;
DataVec=file(2561:2560+NumRows*NumCols*NumFrames);
MaxData=max(DataVec);MinData=min(DataVec);MeanData=mean(DataVec);DataSTD=std(DataVec);

M=zeros(NumCols,NumRows,NumFrames);  % Preallocate matrices
N=M;
M1=zeros(NumCols-2,NumRows-2,NumFrames);
N1=M1;
ActMap=zeros(NumCols,NumRows);
filt=ones(1,NumFrames);
n=filt;
n1=n;

[b,a]=butter(3,[0.05,0.1]); %3rd-order butterworth filter bandpassed near 1Hz
for y = 1:NumRows
    for x = 1:NumCols
        for t=1:NumFrames
            pixel=(y-1)*80+x;
            M(x,y,t)=file((pixel-1)*NumFrames+t+2560);           
       
            mu=mean(M(x,y,:));
            M(x,y,:)=M(x,y,:)-mu;
           
            n(:)=M(x,y,:);
            filt=filter(b,a,n);
            N(x,y,:)=filt;           
           
            %the spatial averaging step (all variables have ones to distingiush
            %from non space-avgd)
            if (y>1) && (y<NumRows) && (x>1) && (x<NumCols)
                x1=x-1;
                y1=y-1;
                M1(x1,y1,:)=M(x,y,:)-(M(x+1,y,:)+M(x-1,y,:)+M(x,y+1,:)+M(x,y-1,:))/4;
                n1(:)=M1(x1,y1,:);
                filt1=filter(b,a,n1);
                N1(x1,y1,:)=filt1;
            end                   
        end        
    end
end

A(INDEX)=M(XCoord,YCoord,INDEX);%raw data at chosen pt.
B(INDEX)=N(XCoord,YCoord,INDEX);%time-filt data at chosen pt.
C(INDEX)=N1(XCoord-1,YCoord-1,INDEX);%time+space filt data at chosen pt.

plot(Times,A);hold on;plot(Times,B,'r','linewidth',2);plot(Times,C,'g','linewidth',2); hold off;
xlabel('Time (s)','FontSize',20)
ylabel('Relative Fluorescence','FontSize',15)
title(['Time Series for ',daFile,' pixel ',num2str(XCoord),',',num2str(YCoord)],'FontSize',15,'Interpreter','none');
%axis([0 NumFrames*FrameInt min(M) max(M)]);
toc