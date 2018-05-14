function [B,A,M,N,NumFrames,NumCols,NumRows,FrameInt,AcquisitionRatio,ibis,ActMap,Times] = daRead4(daFile,XCoord,YCoord)
%this function plots the time series for a single pixel of the redshirt
%camera, both filtered and unfiltered. Pass it the file name and 
%coordinates. (1,1) is top left. It returns three matrices:M, N, ibis
%representing unfiltered data, filtered data, and a binary 3D matrix
%of size M where a 1 denotes an upstroke

%made by Alex Hodge <2008. revised by Bartek Borek in 2008.
tic
data=fopen(daFile);
file=fread(data,'bit16');
NumFrames=file(5);
NumCols=file(385);
NumRows=file(386);
FrameInt=file(389)/1000000;
AcquisitionRatio=file(392);
INDEX = 1:1:NumFrames;
Times = INDEX*FrameInt;
DataVec=file(2561:2560+NumRows*NumCols*NumFrames);
%MaxData=max(DataVec);MinData=min(DataVec);MeanData=mean(DataVec);DataSTD=std(DataVec);
%Diff=MaxData-MinData;


M = zeros(NumCols,NumRows,NumFrames);  % Preallocate matrix
N = M;
ibis=M;
ActMap=zeros(NumCols,NumRows);
filt=ones(1,NumFrames);
n=filt;
[b,a]=butter(3,[0.005,0.2]); %3rd-order butterworth filter bandpassed
%freqz(b,a,128,40) gives you freq response for 40 Hz sampling
for y = 1:NumRows
    for x = 1:NumCols
        for t=1:NumFrames
            pixel=(y-1)*80+x;
            M(x,y,t)=file((pixel-1)*NumFrames+t+2560);           
        end
        mu=mean(M(x,y,:));
        M(x,y,:)=M(x,y,:)-mu; 
        n(:)=M(x,y,:);
        filt=filter(b,a,n);
        N(x,y,:)=filt;
        thresh=(max(filt)+min(filt));
        ActMap(x,y)=(max(filt)-min(filt));
        if (max(N(x,y,:))-min(N(x,y,:)))>600
            for t=1:NumFrames-2
                if N(x,y,t)<thresh
                    if N(x,y,t+1)>=thresh
                        ibis(x,y,t)=1;
                    end
                end               
            end            
        end        
    end
end

A(INDEX)=M(XCoord,YCoord,INDEX);
B(INDEX)=N(XCoord,YCoord,INDEX);
plot(Times,B,'linewidth',2);
hold
plot(Times,A,'r');
xlabel('Time (s)','FontSize',20);
ylabel('Relative Fluorescence','FontSize',15);
title(['Time Series for ',daFile,' pixel ',num2str(XCoord),',',num2str(YCoord)],'FontSize',15,'Interpreter','none');
%axis([0 NumFrames*FrameInt min(M) max(M)]);
toc
