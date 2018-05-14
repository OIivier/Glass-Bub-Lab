function [A,B,C,M1,N1,S1,SF1,NumFrames,NumCols,NumRows,FrameInt,AcquisitionRatio,ActMap,Times] = daRead6(daFile,XCoord,YCoord)
%this function plots the time series for a single pixel of the redshirt
%camera, both filtered and unfiltered. Pass it the file name and 
%coordinates. (1,1) is top left. It returns three matrices:M, N, ibis
%representing unfiltered data, filtered data, and a binary 3D matrix
%of size M where a 1 denotes an upstroke

%made by Alex Hodge <2008. revised by Bartek Borek in 2008.
%daRead5.m reorients so that output matrices are orientied in the same way
%as videos made by IDL cardioplex.
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
%DataVec=file(2561:2560+NumRows*NumCols*NumFrames);
%MaxData=max(DataVec);MinData=min(DataVec);MeanData=mean(DataVec);DataSTD=std(DataVec);
%Diff=MaxData-MinData;


M = zeros(NumCols,NumRows,NumFrames);  % Preallocate matrix
N = M; N1=M; M1=M;
S=zeros(NumCols-2,NumRows-2,NumFrames);
S1=S;SF=S;SF1=S;
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
        ActMap(x,y)=(max(filt)-min(filt));      
    end
end

S=space_avg9(M);

for y = 1:NumRows-2
    for x = 1:NumCols-2
        n(:)=S(x,y,:);
        filt=filter(b,a,n);
        SF(x,y,:)=filt;     
    end
end


for t=1:NumFrames
    M1(:,:,t)=M(:,:,t)';
    N1(:,:,t)=N(:,:,t)';
    S1(:,:,t)=S(:,:,t)';
    SF1(:,:,t)=SF(:,:,t)';
end

A(INDEX)=M1(YCoord,XCoord,INDEX);
B(INDEX)=N1(YCoord,XCoord,INDEX);
C(INDEX)=SF1(YCoord-1,XCoord-1,INDEX);
plot(Times,B,'linewidth',2);
hold
plot(Times,A,'r');
plot(Times,C,'k','linewidth',2);
xlabel('Time (s)','FontSize',15);
ylabel('Relative Fluorescence','FontSize',15);
title(['Time Series for ',daFile,' pixel ',num2str(XCoord),',',num2str(YCoord)],'FontSize',15,'Interpreter','none');
%axis([0 NumFrames*FrameInt min(M) max(M)]);
toc
