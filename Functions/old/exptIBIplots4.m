
%this function plots all IBIs from the different conditions
%for a particular case (the first one I built) see allIBIplot.m
%****remember to name your q so as to remember pixel. ex.: allIBI2060
%(i couldn't make the varname a workspace variable) xD)

allIBIpix1=[spont1_IBIpix1; spont2_IBIpix1; drug1_IBIpix1; drug2_IBIpix1; drug3_IBIpix1; drug4_IBIpix1];
a=length(spont1_IBIpix1);
b=length(spont2_IBIpix1);

beat=1:length(allIBIpix1);

subplot(2,1,1);
plot(beat(1:a+b),allIBIpix1(1:a+b),'b*');
hold;
plot(beat(a+b+1:end),allIBIpix1(a+b+1:end),'r*');
title('IBIs for 13/05/09 dish3 at pix1', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix1=[meanspont1_IBIpix1; meanspont2_IBIpix1; meandrug1_IBIpix1; meandrug2_IBIpix1; meandrug3_IBIpix1; meandrug4_IBIpix1];
allstdpix1=[stdspont1_pix1; stdspont2_pix1; stddrug1_pix1; stddrug2_pix1; stddrug3_pix1; stddrug4_pix1];

subplot(2,1,2);
errorbar([1,2],allmeanpix1(1:2),allstdpix1(1:2),'b*');
hold;
errorbar([3,4,5,6],allmeanpix1(3:6),allstdpix1(3:6),'r*');
xlabel('recording number')
ylabel('mean IBI during recording')
legend('spontaneous', '0.75uM E4031')

hold off

%************pix2
figure

allIBIpix2=[spont1_IBIpix2; spont2_IBIpix2; drug1_IBIpix2; drug2_IBIpix2; drug3_IBIpix2; drug4_IBIpix2];
a=length(spont1_IBIpix2);
b=length(spont2_IBIpix2);

beat=1:length(allIBIpix2);

subplot(2,1,1);
plot(beat(1:a+b),allIBIpix2(1:a+b),'b*');
hold;
plot(beat(a+b+1:end),allIBIpix2(a+b+1:end),'r*');
title('IBIs for 13/05/09 dish3 at pix2', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix2=[meanspont1_IBIpix2; meanspont2_IBIpix2; meandrug1_IBIpix2; meandrug2_IBIpix2; meandrug3_IBIpix2; meandrug4_IBIpix2];
allstdpix2=[stdspont1_pix2; stdspont2_pix2; stddrug1_pix2; stddrug2_pix2; stddrug3_pix2; stddrug4_pix2];

subplot(2,1,2);
errorbar([1,2],allmeanpix2(1:2),allstdpix2(1:2),'b*');
hold;
errorbar([3,4,5,6],allmeanpix2(3:6),allstdpix2(3:6),'r*');
xlabel('recording number')
ylabel('mean IBI during recording')
legend('spontaneous', '0.75uM E4031')

hold off

%************pix3
figure
allIBIpix3=[spont1_IBIpix3; spont2_IBIpix3; drug1_IBIpix3; drug2_IBIpix3; drug3_IBIpix3; drug4_IBIpix3];
a=length(spont1_IBIpix3);
b=length(spont2_IBIpix3);

beat=1:length(allIBIpix3);

subplot(2,1,1);
plot(beat(1:a+b),allIBIpix3(1:a+b),'b*');
hold;
plot(beat(a+b+1:end),allIBIpix3(a+b+1:end),'r*');
title('IBIs for 13/05/09 dish3 at pix3', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix3=[meanspont1_IBIpix3; meanspont2_IBIpix3; meandrug1_IBIpix1; meandrug2_IBIpix1; meandrug3_IBIpix1; meandrug4_IBIpix1];
allstdpix3=[stdspont1_pix3; stdspont2_pix3; stddrug1_pix1; stddrug2_pix1; stddrug3_pix1; stddrug4_pix1];

subplot(2,1,2);
errorbar([1,2],allmeanpix3(1:2),allstdpix3(1:2),'b*');
hold;
errorbar([3,4,5,6],allmeanpix3(3:6),allstdpix3(3:6),'r*');
xlabel('recording number')
ylabel('mean IBI during recording')
legend('spontaneous', '0.75uM E4031')

hold off


%************pix4
figure
allIBIpix4=[spont1_IBIpix4; spont2_IBIpix4; drug1_IBIpix4; drug2_IBIpix4; drug3_IBIpix4; drug4_IBIpix4];
a=length(spont1_IBIpix4);
b=length(spont2_IBIpix4);

beat=1:length(allIBIpix4);

subplot(2,1,1);
plot(beat(1:a+b),allIBIpix4(1:a+b),'b*');
hold;
plot(beat(a+b+1:end),allIBIpix4(a+b+1:end),'r*');
title('IBIs for 13/05/09 dish3 at pix4', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix4=[meanspont1_IBIpix4; meanspont2_IBIpix4; meandrug1_IBIpix4; meandrug2_IBIpix4; meandrug3_IBIpix4; meandrug4_IBIpix4];
allstdpix4=[stdspont1_pix4; stdspont2_pix4; stddrug1_pix4; stddrug2_pix4; stddrug3_pix4; stddrug4_pix4];

subplot(2,1,2);
errorbar([1,2],allmeanpix4(1:2),allstdpix4(1:2),'b*');
hold;
errorbar([3,4,5,6],allmeanpix4(3:6),allstdpix4(3:6),'r*');
xlabel('recording number')
ylabel('mean IBI during recording')
legend('spontaneous', '0.75uM E4031')

hold off