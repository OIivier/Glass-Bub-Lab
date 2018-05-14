% Made by Alex Hodge <2008, revised by Bartek Borek in 2008 and revised 
% again by Lucas Campanari & MinJu You in 3.14.2017. 

function [SF, Lx, Ly] = daRead12(daFile, pixeldx)

% What does it do? 
% 1. Extracts from the da file the data needed to make the video.
% 2. User decides which threshold to use to remove 'obstacles' or 'out of bounds' regions. 
% 3. Resizes matrix based on the pixeldx. Pixeldx = number of columns & rows skipped. 
% 4. Transposes matrix so that output is oriented in the same way as videos made by IDL cardioplex. (1,1) is top left. 

tic

%% Section 1: Extracts data
%daFile = 'MJ4_39.da'; 

data=fopen(daFile);
file=fread(data,'bit16');
NumFrames=file(5);
NumCols=file(385);
NumRows=file(386);
FrameInt=file(389)/1000000;
AcquisitionRatio=file(392);
INDEX = 1:1:NumFrames;
Times = INDEX*FrameInt;

M = zeros(NumCols,NumRows,NumFrames);

for y = 1:NumRows
    for x = 1:NumCols
        for t = 1:NumFrames
            pixel = (y-1)*80+x;
            M(x,y,t) = file((pixel-1)*NumFrames+t+2560);  
        end
        % These next two lines gets the mean of the values at the particular
        % x,y point at all the frames, then subtracts out the mean, to
        % result in a matrix of values that indicate how far the point is 
        % from the mean of that point overtime. 
        mu = mean(M(x,y,:));
        M(x,y,:) = M(x,y,:)-mu;
    end
end



%% Section 2: Resizes matrix M based on pixeldx into newFinal.

[newFinal,Lx,Ly] = ResRed(M,pixeldx);
[nx,ny,NumFrames] = size(newFinal);

%% Optional section: Old method of resizing matrix. Turned on only for comparison. 
% 
% [oriLength,~,~] = size(M);
% dx = pixeldx; % The number of points that will be omitted between two points. 
% 
% nx = ceil(oriLength/(dx+1)); ny = nx; 
% new = zeros(nx, ny, NumFrames);
% newFinal = zeros(nx, ny, NumFrames); 
% 
% % Transferring numbers from matrix M to newFinal 
% for count = 1:NumFrames
%     for k = 1:oriLength
%         oriRowN = M(k,:,count);
%         new(k,:,count) = oriRowN(1:dx+1:end); 
%         % new has desired columns removed. 
%     end
% end
% 
% for count = 1:NumFrames
%     for m = 1:nx
%         newColumnN = new(:,m,count); 
%         newFinal(:,m,count) = newColumnN(1:dx+1:end); 
%         %newFinal has desired columns & rows removed.
%     end
% end
% 
% Lx = 0; 
% Ly = 0; 
%% Section 3: Filtering by Averaging & using Butterworth filter 

% Re-setting variables 
S = zeros(nx,ny,NumFrames); 
SF = S; 

% % Averaging neighbouring pixels 
S = space_avg9_mod(newFinal); 

% Butterworth filter 
ActMap=zeros(nx,ny);
filt=ones(1,NumFrames);
n=filt; % Empty vector of 1s to the length of frames 
[b,a]=butter(3,[0.005,0.4]); %3rd-order butterworth filter bandpassed
%freqz(b,a,128,40) gives you freq response for 40 Hz sampling

for y = 1:nx
   for x = 1:ny
       n(:)=S(x,y,:); % Converts the 3D matrix S(x,y,:) into an array of n that is a vector of pixel x,y at all times. 
       filt=filter(b,a,n); % filters the vector input n with the lowpass function [b,a]
       SF(x,y,:)=filt; % inputs the filtered vector 'filt' into a pixel in SF(x,y)        
   end
end

%% Section 4: Transpose

for t=1:NumFrames
%    SF(:,:,t)=newFinal(:,:,t)';
    SF(:,:,t)=SF(:,:,t)';
end

toc