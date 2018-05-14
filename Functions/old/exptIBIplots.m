function [q]=exptIBIplots(pixel)
%this function plots all IBIs from the different conditions


%set up condition strings
s1='spont1';
c1='cont2';
d1='drug1';
d2='drug2';
%set up combined ibi varname and pixel string
all='allIBI';
pix=num2str(pixel);

%glue them together
s1p=strcat(s1,pix);
c1p=strcat(c1,pix);
d1p=strcat(d1,pix);
d2p=strcat(d2,pix);


q=[spont1_2060; cont2_2060; drug1_2060; drug2_2060];
a=length(spont1_2060);
b=length(cont2_2060);

beat=1:length(allIBI2060);

subplot(2,1,1);
plot(beat(1:a),allIBI2060(1:a),'b*');
hold
plot(beat(a+1:a+b),allIBI2060(a+1:a+b),'g*');
plot(beat(a+b+1:end),allIBI2060(a+b+1:end),'r*');
title('IBIs for 18/04/09 dish4 at pixel (20,60)', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')


allmean2060=[meanspont1_2060; meancont2_2060; meandrug1_2060; meandrug2_2060];
allstd2060=[stdspont1_2060; stdcont2_2060; stddrug1_2060; stddrug2_2060];

subplot(2,1,2);
errorbar(1,allmean2060(1),allstd2060(1),'b*');
hold
errorbar(2,allmean2060(2),allstd2060(2),'g*');
errorbar([3,4],allmean2060(3:4),allstd2060(3:4),'r*');
xlabel('recording')
ylabel('mean IBI during recording')
legend('spontaneous', 'control', '0.75uM E4031')