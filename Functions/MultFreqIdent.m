
%% Finding peak frequency at probes, from bandpass MinJu code
% function rotorboolean = MultFreqIdent(m,Map,fileName)
m = DataFromSimulation; %Comment-swap w/ first line to run non-function
startrec = round(size(m,3)*7/10,0);
m = m(:,:,startrec:end);
%%


nx = size(m,1);
ny = size(m,2);



nProbesX = nx;
nProbesY = ny;

if nProbesX==nx
    dx = 1;
else
    dx = floor(nx/(nProbesX+1));
end

if nProbesY==ny
    dy = 1;
else
       dy = floor(ny/(nProbesY+1));
end

Norm_peakf = zeros(nProbesX,nProbesY);
for xx = 1:nProbesX
    for yy = 1:nProbesY
        pks2= {};
        locs2 = {};
        x = squeeze(m(xx*dx,yy*dy,:));
        [pks,locs] = findpeaks(x);
        for kk = 1:size(pks)
            if pks(kk)>-0.6
                pks2 = [pks2,pks(kk)];
                locs2 = [locs2,locs(kk)];
            end
        end
        z = (length(locs2));
        locs2mat= cell2mat(locs2);
        if (z ~= 0 && z ~= 1)
            Norm_peakf(xx,yy) = 2221/(((locs2mat(z))-locs2mat(1))/z);
        elseif z == 1
            Norm_peakf(xx,yy) = 1/3;
        end
    end
end
    
        

%% Deciding if rotor or not
fig=figure('Position', [100, 100, 500, 800]);

subplot(2,1,1);
Heatmap3(Norm_peakf,jet,0.0,4.0);

subplot(2,1,2);
histogram(Norm_peakf);%, 'BinWidth', 3e-5);
xlabel('Frequency (Hz)')
ylabel('Count')
set(gca, 'YScale', 'log')
% 
% subplot(2,2,3);
% Heatmap3(int16(Map),hsv,-1,1);
% xlim([0.5 5]);
% ylim([0 20]);
% 
% 
% print(fig,fileName,'-dpng')
% print(fig,fileName,'-depsc')

% rotorboolean = 1;
% error = 0.15;
% plot(squeeze(m(4*dx,4*dy,:)));
% disp(Norm_peakf);
% averagef = 0;
% counter= 0;
% for x = 1:nProbesX
%     for y = 1:nProbesY
%         if Norm_peakf(x,y) ~= 0
%             averagef = averagef+Norm_peakf(x,y);
%             counter=counter+1;
%         end
%     end
% end
% averagef = averagef/counter;
% problems = 0;
% for xx = 1:nProbesX
%     for yy = 1:nProbesY
%         if Norm_peakf(xx,yy) ~= 0
%             if abs((Norm_peakf(xx,yy)-averagef))/averagef > error
%                 problems = problems+1;
%                 fprintf('Problem at pixel (%d,%d)\n', xx,yy);
%             end
%         end
%     end
% end
% fprintf('problems = %d\n', problems);
% % fprintf('averagef = '+averagef);
% if isequal(Norm_peakf,zeros(nProbesX,nProbesY))||problems>3
%     rotorboolean = 0;
% end
% RotorPresence = {'No Rotor', 'Yes Rotor'};
% fprintf('%s\n', RotorPresence{1+rotorboolean});
% std(std(Norm_peakf))
% plot(squeeze(m(dx*5,dy*6,:)));
Norm_peakf