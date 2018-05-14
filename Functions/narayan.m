%%
% Change this to the path in your computer where this script is
cd('/home/vincent/Documents/McGill/PHYS 459 - Research Thesis/Data')

fileID = fopen('Case 1 QRS Subtracted Basket Data.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);

%
N = 4000; % number of frames
Asz = size(A,1);
P = reshape(A,N,Asz/N);

sP = std(P,0,1);
mP = mean(P,1);

% Standard deviation of each probe
figure
plot(sP,'-ob','linewidth', 1.5)
ylabel 'Std'
xlabel 'Probe number'
axis tight

% All probes
figure
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01]);

for ii = 1 : 64  
  axes(ha(ii)); 
  plot(P(:,ii), 'linewidth', 1.1);

  text(0.5,0.9,num2str(ii),'units','normalized','fontweight','bold')  
  set(gca,'xtick',[])
  set(gca,'xticklabel',[],'fontsize',8)
%   set(gca,'ytick',[])
%   set(gca,'yticklabel',[])
axis tight

end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
