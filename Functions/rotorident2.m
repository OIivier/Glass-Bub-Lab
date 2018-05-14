
%% Finding peak frequency at probes, from bandpass MinJu code
function rotorboolean = rotorident2(m)
temp = m;
nx = size(temp,1);
ny = size(temp,2);

nProbesX = 40;
nProbesY = 40;

dx = floor(nx/(nProbesX+1));
dy = floor(ny/(nProbesY+1));
startrec = round(size(temp,3)*9/10,0);
Norm_peakf = zeros(nProbesX,nProbesY);
for ii = 1:nProbesX
    for jj = 1:nProbesY
        pks2= {};
        locs2 = {};
        x = squeeze(m(ii*dx,jj*dy,startrec:end));
        [pks,locs] = findpeaks(x);
        for kk = 1:size(pks)
            if pks(kk)>0
                pks2 = [pks2,pks(kk)];
                locs2 = [locs2,locs(kk)];
            end
        end
        z = (length(locs2));
        locs2mat= cell2mat(locs2);
        if z ~= 0
            Norm_peakf(ii,jj) = 2221/(((locs2mat(z))-locs2mat(1))/z);
        end
    end
end
    
        

%% Deciding if rotor or not

rotorboolean = 1;
error = 0.15;
% plot(squeeze(m(4*dx,4*dy,:)));
disp(Norm_peakf);
averagef = 0;
counter= 0;
for x = 1:nProbesX
    for y = 1:nProbesY
        if Norm_peakf(x,y) ~= 0
            averagef = averagef+Norm_peakf(x,y);
            counter=counter+1;
        end
    end
end
averagef = averagef/counter;
problems = 0;
for j = 1:nProbesX
    for k = 1:nProbesY
        if Norm_peakf(j,k) ~= 0
            if abs((Norm_peakf(j,k)-averagef))/averagef > error
                problems = problems+1;
            end
        end
    end
end
fprintf('problems = %d \n', problems);
% fprintf('averagef = '+averagef);
if isequal(Norm_peakf,zeros(nProbesX,nProbesY))||problems>3
    rotorboolean = 0;
end
RotorPresence = {'No Rotor', 'Yes Rotor'};
fprintf('%s\n', RotorPresence{1+rotorboolean});
std(std(Norm_peakf))
% plot(squeeze(m(dx*5,dy*6,:)));
plot(squeeze(m(dx*23,dy*30,:)));
hold on 
plot(squeeze(m(dx*30,dy*23,:)));