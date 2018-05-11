function [phaseMap] = phaseMap_1_old(SF2)

% This code creates phasemap of the data. Then the angles are corrected so
% that the leading edge has the lowest value (-180~) and trailing edge has
% the highest value (~+180). For example, if you were to color it with 5
% colors, -180~-108 = red, and +108~180 = blue. 

tic
%%
% MAKE SURE YOU PUT IN THE RIGHT FILE (ie. SF1, SF2, SF??) 
temp = SF2;

[row, column, frame] = size(temp);
phaseMap = zeros(row, column, frame); 

iniT = 1; 
finT = frame;
delay = 5;  %Delay is normally given a value that is 1/4 of the Period. For the time being, I've given it an arbitrary value.

% Create a matrix with values corresponding to the phases at each pixel  
for n = 1:row 
    for m = 1:column 
        current1Raw = squeeze(temp(n,m,iniT:finT-delay)); %g(t)
        current2Raw = squeeze(temp(n,m,iniT+delay:finT)); %g(t+1)
        mean1 = 0; %mean(current1Raw);
        mean2 = 0; %mean(current2Raw);
        
        for k = 1:frame-delay
            rawAngle = atan2d((current2Raw(k)-mean2),(current1Raw(k)-mean1)); % This line calculates the angle.
            
            % The following lines 'corrects' the values of the angles so
            % that they run from -180 to +180, red to blue. These lines are
            % necessary for the correct visual representation using matlab.
            NegAngle = rawAngle*-1;
            if NegAngle < -90
               angle = NegAngle + 270;
            else 
               angle = NegAngle - 90; 
            end
            phaseMap(n,m,k) = angle;
        end
    end
end

toc