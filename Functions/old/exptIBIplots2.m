
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
title('IBIs for 20/03/09 dish2 at pixel (27,31)', 'Fontsize',14)
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
legend('spontaneous', 'control', '0.75uM E4031')

hold off

%************pix2
figure

allIBIpix2=[spont1_IBIpix2; cont1_IBIpix2; drug1_IBIpix2; drug2_IBIpix2];
a=length(spont1_IBIpix2);
b=length(cont1_IBIpix2);

beat=1:length(allIBIpix2);

subplot(2,1,1);
plot(beat(1:a),allIBIpix2(1:a),'b*');
hold
plot(beat(a+1:a+b),allIBIpix2(a+1:a+b),'g*');
plot(beat(a+b+1:end),allIBIpix2(a+b+1:end),'r*');
title('IBIs for 20/03/09 dish2 at pixel (20,65)', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix2=[meanspont1_IBIpix2; meancont1_IBIpix2; meandrug1_IBIpix2; meandrug2_IBIpix2];
allstdpix2=[stdspont1_pix2; stdcont1_pix2; stddrug1_pix2; stddrug2_pix2];

subplot(2,1,2);
errorbar(1,allmeanpix2(1),allstdpix2(1),'b*');
hold
errorbar(2,allmeanpix2(2),allstdpix2(2),'g*');
errorbar([3,4],allmeanpix2(3:4),allstdpix2(3:4),'r*');
xlabel('recording')
ylabel('mean IBI during recording')
legend('spontaneous', 'control', '0.75uM E4031')

hold off

%************pix3
figure

allIBIpix3=[spont1_IBIpix3; cont1_IBIpix3; drug1_IBIpix3; drug2_IBIpix3];
a=length(spont1_IBIpix3);
b=length(cont1_IBIpix3);

beat=1:length(allIBIpix3);

subplot(2,1,1);
plot(beat(1:a),allIBIpix3(1:a),'b*');
hold
plot(beat(a+1:a+b),allIBIpix3(a+1:a+b),'g*');
plot(beat(a+b+1:end),allIBIpix3(a+b+1:end),'r*');
title('IBIs for 20/03/09 dish2 at pixel (32,51)', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix3=[meanspont1_IBIpix3; meancont1_IBIpix3; meandrug1_IBIpix3; meandrug2_IBIpix3];
allstdpix3=[stdspont1_pix3; stdcont1_pix3; stddrug1_pix3; stddrug2_pix3];

subplot(2,1,2);
errorbar(1,allmeanpix3(1),allstdpix3(1),'b*');
hold
errorbar(2,allmeanpix3(2),allstdpix3(2),'g*');
errorbar([3,4],allmeanpix3(3:4),allstdpix3(3:4),'r*');
xlabel('recording')
ylabel('mean IBI during recording')
legend('spontaneous', 'control', '0.75uM E4031')

hold off

%************pix4
figure

allIBIpix4=[spont1_IBIpix4; cont1_IBIpix4; drug1_IBIpix4; drug2_IBIpix4];
a=length(spont1_IBIpix4);
b=length(cont1_IBIpix4);

beat=1:length(allIBIpix4);

subplot(2,1,1);
plot(beat(1:a),allIBIpix4(1:a),'b*');
hold
plot(beat(a+1:a+b),allIBIpix4(a+1:a+b),'g*');
plot(beat(a+b+1:end),allIBIpix4(a+b+1:end),'r*');
title('IBIs for 20/03/09 dish2 at pixel (63,66)', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix4=[meanspont1_IBIpix4; meancont1_IBIpix4; meandrug1_IBIpix4; meandrug2_IBIpix4];
allstdpix4=[stdspont1_pix4; stdcont1_pix4; stddrug1_pix4; stddrug2_pix4];

subplot(2,1,2);
errorbar(1,allmeanpix4(1),allstdpix4(1),'b*');
hold
errorbar(2,allmeanpix4(2),allstdpix4(2),'g*');
errorbar([3,4],allmeanpix4(3:4),allstdpix4(3:4),'r*');
xlabel('recording')
ylabel('mean IBI during recording')
legend('spontaneous', 'control', '0.75uM E4031')
