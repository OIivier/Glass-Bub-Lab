Data = OpenNarayan(6);
MaximumArray = zeros(64,1);
DataAutocor = Data;
Scale = 2;
AutocorBool = 1;
if AutocorBool==1
    for i=1:64  
        DataAutocor(:,i) = AutocorHM(Data(:,i))';
        plot(Data(:,1));
        plot(DataAutocor(1:end-100));
        DataFT = abs(fft(DataAutocor(:,1)));
        DataFT = DataFT(1:floor(size(DataAutocor,1)/2));
        [M,I] = max(DataFT);
        MaximumArray(i) = I;
    end
end
    
for i=1:8
    for j=1:8
        subplot(8,8,(i-1)*8+j);
        plot(squeeze(DataAutocor(:,(i-1)*8+j)));
        title(['(' num2str(i) ',' num2str(j) ')'])
        xlabel('Time (millisecond)')
        ylabel('Voltage (V)')
%         ylim([-Scale Scale])
        ax = gca;
        ax.FontSize = 6;
    end
end



% subplot(2,2,1);
% plot(squeeze(DataAutocor(:,2)));
% title(['Probe (1,2)'])
% xlabel('Time (millisecond)')
% ylabel('Voltage (mV)')
% ylim([-Scale Scale])
% ax = gca;
% ax.FontSize = 12;
% 
% subplot(2,2,2);
% plot(squeeze(DataAutocor(:,5)));
% title(['Probe (1,5)'])
% xlabel('Time (millisecond)')
% ylabel('Voltage (mV)')
% ylim([-Scale Scale])
% ax = gca;
% ax.FontSize = 12;
% 
% subplot(2,2,3);
% plot(squeeze(DataAutocor(:,54)));
% title(['Probe (7,6)'])
% xlabel('Time (millisecond)')
% ylabel('Voltage (mV)')
% ylim([-Scale Scale])
% ax = gca;
% ax.FontSize = 12;
% 
% subplot(2,2,4);
% plot(squeeze(DataAutocor(:,59)));
% title(['Probe (8,4)'])
% xlabel('Time (millisecond)')
% ylabel('Voltage (mV)')
% ylim([-Scale Scale])
% ax = gca;
% ax.FontSize = 12;

set(gcf, 'Position', [0, 0, 1200, 1200])
















% count=0;
% for i=1:64
%     if MaximumArray(i)>5
%        count=count+1;
%     end
% end

% j=1;
% MaximumArrayGood = zeros(count,1);
% for i=1:64
%     if MaximumArray(i)>5
%         MaximumArrayGood(j)=MaximumArray(i)/4000;
%         j=j+1;
%     end
% end

% histogram(MaximumArrayGood, 'BinWidth', 3e-5);
% xlim([0 8e-3]);
% ylim([0 20]);
% xlabel('Frequency (frame^{-1})')
% ylabel('Count')