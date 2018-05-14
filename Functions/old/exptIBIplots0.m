
%this function plots all IBIs from the different conditions
%for a particular case (the first one I built) see allIBIplot.m
%****remember to name your q so as to remember pixel. ex.: allIBI2060
%(i couldn't make the varname a workspace variable) xD)

allIBIpix1=[spont1_IBIpix1; cont1_IBIpix1; drug1_IBIpix1; drug2_IBIpix1];
a=length(spont1_IBIpix1);
b=length(cont1_IBIpix1);

beat=1:length(allIBIpix1);

subplot(2,1,1);
plot(beat(1:a),allIBIpix1(1:a),'b*');
hold
plot(beat(a+1:a+b),allIBIpix1(a+1:a+b),'g*');
plot(beat(a+b+1:end),allIBIpix1(a+b+1:end),'r*');
title('IBIs for 08/04/09 dish3 at pix1', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix1=[meanspont1_IBIpix1; meancont1_IBIpix1; meandrug1_IBIpix1; meandrug2_IBIpix1];
allstdpix1=[stdspont1_pix1; stdcont1_pix1; stddrug1_pix1; stddrug2_pix1];

subplot(2,1,2);
errorbar(1,allmeanpix1(1),allstdpix1(1),'b*');
hold
errorbar(2,allmeanpix1(2),allstdpix1(2),'g*');
errorbar([3,4],allmeanpix1(3:4),allstdpix1(3:4),'r*');
xlabel('recording')
ylabel('mean IBI during recording')
legend('spontaneous', 'control', '1.5uM E4031')

hold off
