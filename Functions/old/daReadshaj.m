function [SF1,ActMap,Times] = daRead7(daFile)
%this function plots the time series for a single pixel of the redshirt
%camera, both filtered and unfiltered. Pass it the file name and 
%coordinates. (1,1) is top left. It returns the matrix SF which have each
%pixel space averaged and time filtered. 

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

% Preallocate space for matrices
M = zeros(NumCols,NumRows,NumFrames); 
S=zeros(NumCols-2,NumRows-2,NumFrames);
SF=S;SF1=S;
ActMap=zeros(NumCols,NumRows);
filt=ones(1,NumFrames);
n=filt;
[b,a]=butter(3,[0.005,0.4]); %3rd-order butterworth filter bandpassed
%freqz(b,a,128,40) gives you freq response for 40 Hz sampling

%shaj
NumFrames2=NumFrames/10;
for i=1:10

for y = 1:NumRows
    for x = 1:NumCols
        for t=1:NumFrames2
            pixel=(y-1)*80+x;
            M(x,y,t)=file((pixel-1)*NumFrames+(i-1)*NumFrames2+t+2560);  
        end
        mu=mean(M(x,y,:));
        M(x,y,:)=M(x,y,:)-mu;             
    end
end

%now we average each cell with 8 of it's nearest neighbours.**** remember
%S(1,1) represents M(2,2)!
S=space_avg9(M);

for y = 1:NumRows-2
    for x = 1:NumCols-2
        n(:)=S(x,y,:);
        filt=filter(b,a,n);
        SF(x,y,:)=filt;  
        ActMap(x,y)=(max(filt)-min(filt));     
    end
end


for t=1:NumFrames2
    SF1(:,:,t)=SF(:,:,t);
colormap('summer')
cmin=min(SF1(20,41,:));
cmax=max(SF1(20,41,:));
pcolor(SF1(:,:,t));shading interp;colorbar;
caxis([cmin,cmax])
k=(i-1)*NumFrames2+t
	if (k<10)
	print('-dpng', ['plot_000', int2str(k)]);
	elseif (k<100)
	print('-dpng', ['plot_00', int2str(k)]);
	elseif (k<1000)
	print('-dpng', ['plot_0', int2str(k)]);
	else 
	print('-dpng', ['plot_', int2str(k)]);
	end

end

end

toc
