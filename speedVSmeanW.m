wLStart = 0.2;
wLEnd = 0.8;
nPoints = 20;
wLStepSize = (wLEnd-wLStart)/nPoints;
% meanWArray = zeros(nPoints,1);
% speedArray = zeros(nPoints,1);
% lengthArray = zeros(nPoints,1);
endDistanceArray = zeros(nPoints,1);
initDistanceArray = zeros(nPoints,1);
cycleLengthArray = zeros(nPoints,1);
apdAArray=zeros(nPoints,1);
apdBArray=zeros(nPoints,1);
apdCArray=zeros(nPoints,1);
RealLength=1;
dx = 0.005;
RealDuration=6;
dt = 0.01;
DataPointArray = zeros(nPoints,round(RealDuration/(dt^2)));
startLength=0.05;
endLength=1.2;

for i=1:nPoints
%     [cycleLengthArray(i) apdAArray(i) apdCArray(i) apdBArray(i) DataPointArray(i,:,:)] = Sim1D(wLStart+(i-1)*wLStepSize,0,10000,RealLength, dx, RealDuration, dt, 90,0.1);
%     disp(i);
    endDistanceArray(i) = SimRing(0.54, 0, startLength+(i-1)*(endLength-startLength)/nPoints, 1, dx,2,dt,20);
    initDistanceArray(i) = startLength+(i-1)*(endLength-startLength)/nPoints;
end
    
fprintf('\n Done\n');

%%

% hold off


figure;



% plot(cycleLengthArray-apdAArray,apdAArray,'LineWidth',2);
% hold on
% plot(cycleLengthArray-apdBArray,apdBArray,'LineWidth',2);
% plot(cycleLengthArray-apdCArray,apdCArray,'LineWidth',2);
% xlabel('Recovery Period');
% ylabel('APD');

plot(initDistanceArray,2*dx*round(endDistanceArray/(2*dx)),'LineWidth',2);
xlabel('Initial distance (cm)');
ylabel('End distance (cm)');
cat(2,initDistanceArray,endDistanceArray)
% 
% ylabel('Final distance (px)');
% subplot(3,1,1);
% plot(meanWArray,speedArray,'LineWidth',2);
% xlabel('Mean w');
% xlim([min(min(meanWArray)),max(max(meanWArray))]);
% ylabel('Speed');
% subplot(3,1,2);
% plot(lengthArray,speedArray,'LineWidth',2);
% xlabel('Length');
% ylabel('Speed');
% subplot(3,1,3);
% plot(lengthArray,meanWArray,'LineWidth',2);
% xlabel('Length');
% ylabel('Mean w');

%%

figure('pos',[200 200 1000 1000]);

ax = gca;
ax.NextPlot = 'replaceChildren';

cd('/home/vincent/Documents/Glass-Bub-Lab/FiguresAndGifs');


nFrames = size(DataPointArray,2);
F(nFrames) = struct('cdata',[],'colormap',[]);

firstPeak = 10;
for j = 1:round(nFrames/1)
    [peaks locs] = findpeaks(squeeze(DataPointArray(j,:)));
    plot(DataPointArray(j,500*floor(locs(firstPeak)/500):500*ceil(locs(firstPeak+2)/500)),'LineWidth',2);
%     xlim([locs(5) locs(7)]);
%     s = surf(DataSkipFrames(:,:,j)');
%     s.EdgeColor = 'none';
%     caxis([-1.95 1.7]);
%     view([0 90]);




    drawnow

    F(j) = getframe;

    im = frame2im(F(j)); 
    [imind,cm] = rgb2ind(im,128); 

    %Write to the GIF File 
%     if j == 1 
%         imwrite(imind,cm,file.name,'gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,file.name,'gif','WriteMode','append'); 
%     end 
end