figure('pos',[100 500 1000 1000]);

ax = gca;
ax.NextPlot = 'replaceChildren';

cd('/home/vincent/Documents/Research Thesis/FiguresAndGifs');

% Data = Data(4:Length+3,4,1:40:end);
Diameter = ceil(size(Data,1)/pi);
RingMap = zeros(Diameter+2,Diameter+2,size(Data,3));

for i=1:size(Data,1)

    Angle = i*2*pi/size(Data,1);
    X = round((Diameter+2)/2 + (Diameter/2)*cos(Angle)); 
    Y = round((Diameter+2)/2 + (Diameter/2)*sin(Angle));
    RingMap(X,Y,:)=Data(i,:);
    
end
    

nFrames = size(Data,3);
F(nFrames) = struct('cdata',[],'colormap',[]);


for j = 1:round(nFrames/1)
    s = surf(RingMap(:,:,j)');
    s.EdgeColor = 'none';
    caxis([-1.95 1.7]);
    view([0 90]);
        
    
       
    
    drawnow
    
    F(j) = getframe;
      
    im = frame2im(F(j)); 
    [imind,cm] = rgb2ind(im,128); 
    
    %Write to the GIF File 
    if j == 1 
        imwrite(imind,cm,file.name,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,file.name,'gif','WriteMode','append'); 
    end 
end