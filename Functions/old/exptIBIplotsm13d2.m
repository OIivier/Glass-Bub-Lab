%this function plots all IBIs from the different conditions
%for a particular case (the first one I built) see allIBIplot.m
%****remember to name your q so as to remember pixel. ex.: allIBI2060
%(i couldn't make the varname a workspace variable) xD)
close all;
%% pix1 %%
figure(1)
allIBIpix1=[spont_noperf_IBIpix1; spont1_perf_IBIpix1; spont2_perf_IBIpix1; spont3_perf_IBIpix1; drug1_IBIpix1; drug2_IBIpix1; drug3_IBIpix1; drug4_IBIpix1;];

a=length(spont_noperf_IBIpix1);
b=length(spont1_perf_IBIpix1);
c=length(spont2_perf_IBIpix1)+length(spont3_perf_IBIpix1)+length(drug1_IBIpix1);


beat=1:length(allIBIpix1);

subplot(4,1,1);
plot(beat(1:a),allIBIpix1(1:a),'g*');
hold;
plot(beat(a+1:a+b),allIBIpix1(a+1:a+b),'b*');
plot(beat(a+b+1:a+b+c),allIBIpix1(a+b+1:a+b+c),'g*');
plot(beat(a+b+c+1:end),allIBIpix1(a+b+c+1:end),'b*')
title('IBIs for 13/05/09 dish2 at (40,60)', 'Fontsize',30)
ylabel('IBI [s]','FontSize',16)
xlabel('Beat Number Analyzed','FontSize',16)
ylim([0 2.5]);

allmeanpix1=[meanspont_noperf_IBIpix1; meanspont1_perf_IBIpix1; meanspont2_perf_IBIpix1; meanspont3_perf_IBIpix1; meandrug1_IBIpix1; meandrug2_IBIpix1; meandrug3_IBIpix1; meandrug4_IBIpix1;];
allstdpix1=[stdspont_noperf_pix1; stdspont1_perf_pix1; stdspont2_perf_pix1; stdspont3_perf_pix1; stddrug1_pix1; stddrug2_pix1; stddrug3_pix1; stddrug4_pix1;];

subplot(4,1,2);
errorbar([1],allmeanpix1(1),allstdpix1(1),'g*');
hold;
errorbar([2],allmeanpix1(2),allstdpix1(2),'b*');
errorbar([3:5],allmeanpix1(3:5),allstdpix1(3:5),'g*');
errorbar([6:8],allmeanpix1(6:8),allstdpix1(6:8),'b*');
xlim([1 8.5]);
grid on;
%xlabel('Recording number','FontSize',16);
ylabel('Mean IBI During Recording','FontSize',16);
legend('Spontaneous (No Perfusion) (0:00)','Spontaneous (Perfusion (1.15 mL/min)) (0:03)', 'Spontaneous (No Perfusion) (0:13-0:40)', 'Spontaneous (Perfusion (0.5 mL/min,0.75mL/min,1.0mL/min) (0:51-1:08))','FontSize','16');

subplot(4,1,3);
recnum=1:8;
Therm=[35 35 36.5 36.5 34.5 34.5 34 33.5];
RealTemp=[36 36 38 38 35.5 35.5 35 34.5];
msize=15;
plot(recnum(1),Therm(1),'g*','MarkerSize',msize); ylim([33 39]);
grid on;
xlim([1 8.5])
hold;
plot(recnum(1),RealTemp(1),'go','MarkerSize',msize);
plot(recnum(2),Therm(2),'b*','MarkerSize',msize);
plot(recnum(2),RealTemp(2),'bo','MarkerSize',msize);
plot(recnum(3:5),Therm(3:5),'g*','MarkerSize',msize);
plot(recnum(3:5),RealTemp(3:5),'go','MarkerSize',msize);
plot(recnum(6:8),Therm(6:8),'b*','MarkerSize',msize);
plot(recnum(6:8),RealTemp(6:8),'bo','MarkerSize',msize);
ylabel('Temperature','FontSize',16);

subplot(4,1,4);
perfrate=[0 1.15 0 0 0 1.15 1.725 2.3];
plot(recnum(1),perfrate(1),'g*');
hold on;
plot(recnum(2),perfrate(2),'b*');
plot(recnum(3:5),perfrate(3:5),'g*');
plot(recnum(6:8),perfrate(6:8),'b*');
grid on;
xlim([1 8.5])
xlabel('Recording Number','FontSize',16); ylabel('Perfusion Rate (mL/min)','FontSize',16);

hold off

%% pix2 %%
figure(2)
allIBIpix2=[spont_noperf_IBIpix2; spont1_perf_IBIpix2; spont2_perf_IBIpix2; spont3_perf_IBIpix2; drug1_IBIpix2; drug2_IBIpix2; drug3_IBIpix2; drug4_IBIpix2;];

a=length(spont_noperf_IBIpix2);
b=length(spont1_perf_IBIpix2);
c=length(spont2_perf_IBIpix2)+length(spont3_perf_IBIpix2)+length(drug1_IBIpix2);


beat=1:length(allIBIpix2);

subplot(4,1,1);
plot(beat(1:a),allIBIpix2(1:a),'g*');
hold;
plot(beat(a+1:a+b),allIBIpix2(a+1:a+b),'b*');
plot(beat(a+b+1:a+b+c),allIBIpix2(a+b+1:a+b+c),'g*');
plot(beat(a+b+c+1:end),allIBIpix2(a+b+c+1:end),'b*')
title('IBIs for 13/05/09 dish2 at (40,40)', 'Fontsize',30)
ylabel('IBI [s]', 'Fontsize',16)
xlabel('Beat Number Analyzed', 'Fontsize',16)
ylim([0 2.5]);

allmeanpix2=[meanspont_noperf_IBIpix2; meanspont1_perf_IBIpix2; meanspont2_perf_IBIpix2; meanspont3_perf_IBIpix2; meandrug1_IBIpix2; meandrug2_IBIpix2; meandrug3_IBIpix2; meandrug4_IBIpix2;];
allstdpix2=[stdspont_noperf_pix2; stdspont1_perf_pix2; stdspont2_perf_pix2; stdspont3_perf_pix2; stddrug1_pix2; stddrug2_pix2; stddrug3_pix2; stddrug4_pix2;];

subplot(4,1,2);
errorbar([1],allmeanpix2(1),allstdpix2(1),'g*');
hold;
errorbar([2],allmeanpix2(2),allstdpix2(2),'b*');
errorbar([3:5],allmeanpix2(3:5),allstdpix2(3:5),'g*');
errorbar([6:8],allmeanpix2(6:8),allstdpix2(6:8),'b*');
xlim([1 8.5])
ylabel('Mean IBI During Recording', 'Fontsize',16);
legend('Spontaneous (No Perfusion) (0:00)','Spontaneous (Perfusion (0.5 mL/min)) (0:03)', 'Spontaneous (No Perfusion) (0:13-0:40)', 'Spontaneous (Perfusion (0.5 mL/min,0.75mL/min,1.0mL/min) (0:51-1:08))');

subplot(4,1,3);
recnum=1:8;
Therm=[35 35 36.5 36.5 34.5 34.5 34 33.5];
RealTemp=[36 36 38 38 35.5 35.5 35 34.5];
msize=15;
plot(recnum(1),Therm(1),'g*','MarkerSize',msize); ylim([33 39]);
grid on;
xlim([1 8.5])
hold;
plot(recnum(1),RealTemp(1),'go','MarkerSize',msize);
plot(recnum(2),Therm(2),'b*','MarkerSize',msize);
plot(recnum(2),RealTemp(2),'bo','MarkerSize',msize);
plot(recnum(3:5),Therm(3:5),'g*','MarkerSize',msize);
plot(recnum(3:5),RealTemp(3:5),'go','MarkerSize',msize);
plot(recnum(6:8),Therm(6:8),'b*','MarkerSize',msize);
plot(recnum(6:8),RealTemp(6:8),'bo','MarkerSize',msize);
ylabel('Temperature','FontSize',16);

subplot(4,1,4);
perfrate=[0 1.15 0 0 0 1.15 1.725 2.3];
plot(recnum(1),perfrate(1),'g*','Markersize',msize);
hold on;
plot(recnum(2),perfrate(2),'b*','Markersize',msize);
plot(recnum(3:5),perfrate(3:5),'g*','Markersize',msize);
plot(recnum(6:8),perfrate(6:8),'b*','Markersize',msize);
grid on;
xlim([1 8.5])
xlabel('Recording Number','FontSize',16); ylabel('Perfusion Rate (mL/min)','FontSize',16);

hold off

%% pix3 %%
figure(3)
allIBIpix3=[spont_noperf_IBIpix3; spont1_perf_IBIpix3; spont2_perf_IBIpix3; spont3_perf_IBIpix3; drug1_IBIpix3; drug2_IBIpix3; drug3_IBIpix3; drug4_IBIpix3;];

a=length(spont_noperf_IBIpix3);
b=length(spont1_perf_IBIpix3);
c=length(spont2_perf_IBIpix3)+length(spont3_perf_IBIpix3)+length(drug1_IBIpix3);


beat=1:length(allIBIpix3);

subplot(4,1,1);
plot(beat(1:a),allIBIpix3(1:a),'g*');
hold;
plot(beat(a+1:a+b),allIBIpix3(a+1:a+b),'b*');
plot(beat(a+b+1:a+b+c),allIBIpix3(a+b+1:a+b+c),'g*');
plot(beat(a+b+c+1:end),allIBIpix3(a+b+c+1:end),'b*')
title('IBIs for 13/05/09 dish2 at (30,50)', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix3=[meanspont_noperf_IBIpix3; meanspont1_perf_IBIpix3; meanspont2_perf_IBIpix3; meanspont3_perf_IBIpix3; meandrug1_IBIpix3; meandrug2_IBIpix3; meandrug3_IBIpix3; meandrug4_IBIpix3;];
allstdpix3=[stdspont_noperf_pix3; stdspont1_perf_pix3; stdspont2_perf_pix3; stdspont3_perf_pix3; stddrug1_pix3; stddrug2_pix3; stddrug3_pix3; stddrug4_pix3;];

subplot(4,1,2);
errorbar([1],allmeanpix3(1),allstdpix3(1),'g*');
hold;
errorbar([2],allmeanpix3(2),allstdpix3(2),'b*');
errorbar([3:5],allmeanpix3(3:5),allstdpix3(3:5),'g*');
errorbar([6:8],allmeanpix3(6:8),allstdpix3(6:8),'b*');
xlabel('recording number');
ylabel('mean IBI during recording');
legend('Spontaneous (No Perfusion) (0:00)','Spontaneous (Perfusion (0.5 mL/min)) (0:03)', 'Spontaneous (No Perfusion) (0:13-0:40)', 'Spontaneous (Perfusion (0.5 mL/min,0.75mL/min,1.0mL/min) (0:51-1:08))');

subplot(4,1,3);
recnum=1:8;
Therm=[35 35 36.5 36.5 34.5 34.5 34 33.5];
RealTemp=[36 36 38 38 35.5 35.5 35 34.5];
msize=15;
plot(recnum(1),Therm(1),'g*','MarkerSize',msize); ylim([33 39]);
grid on;
xlim([1 8.5])
hold;
plot(recnum(1),RealTemp(1),'go','MarkerSize',msize);
plot(recnum(2),Therm(2),'b*','MarkerSize',msize);
plot(recnum(2),RealTemp(2),'bo','MarkerSize',msize);
plot(recnum(3:5),Therm(3:5),'g*','MarkerSize',msize);
plot(recnum(3:5),RealTemp(3:5),'go','MarkerSize',msize);
plot(recnum(6:8),Therm(6:8),'b*','MarkerSize',msize);
plot(recnum(6:8),RealTemp(6:8),'bo','MarkerSize',msize);
ylabel('Temperature','FontSize',16);

subplot(4,1,4);
perfrate=[0 1.15 0 0 0 1.15 1.725 2.3];
plot(recnum(1),perfrate(1),'g*');
hold on;
plot(recnum(2),perfrate(2),'b*');
plot(recnum(3:5),perfrate(3:5),'g*');
plot(recnum(6:8),perfrate(6:8),'b*');
grid on;
xlim([1 8.5])
xlabel('Recording Number','FontSize',16); ylabel('Perfusion Rate (mL/min)','FontSize',16);

hold off

%% pix4 %% 
figure(4)
allIBIpix4=[spont_noperf_IBIpix4; spont1_perf_IBIpix4; spont2_perf_IBIpix4; spont3_perf_IBIpix4; drug1_IBIpix4; drug2_IBIpix4; drug3_IBIpix4; drug4_IBIpix4;];

a=length(spont_noperf_IBIpix4);
b=length(spont1_perf_IBIpix4);
c=length(spont2_perf_IBIpix4)+length(spont3_perf_IBIpix4)+length(drug1_IBIpix4);


beat=1:length(allIBIpix4);

subplot(4,1,1);
plot(beat(1:a),allIBIpix4(1:a),'g*');
hold;
plot(beat(a+1:a+b),allIBIpix4(a+1:a+b),'b*');
plot(beat(a+b+1:a+b+c),allIBIpix4(a+b+1:a+b+c),'g*');
plot(beat(a+b+c+1:end),allIBIpix4(a+b+c+1:end),'b*')
title('IBIs for 13/05/09 dish2 at (50,40)', 'Fontsize',14)
ylabel('IBI [s]')
xlabel('beat number analyzed')
ylim([0 2.5]);

allmeanpix4=[meanspont_noperf_IBIpix4; meanspont1_perf_IBIpix4; meanspont2_perf_IBIpix4; meanspont3_perf_IBIpix4; meandrug1_IBIpix2; meandrug2_IBIpix4; meandrug3_IBIpix4; meandrug4_IBIpix4;];
allstdpix4=[stdspont_noperf_pix4; stdspont1_perf_pix4; stdspont2_perf_pix4; stdspont3_perf_pix4; stddrug1_pix4; stddrug2_pix4; stddrug3_pix4; stddrug4_pix4;];

subplot(4,1,2);
errorbar([1],allmeanpix4(1),allstdpix4(1),'g*');
hold;
errorbar([2],allmeanpix4(2),allstdpix4(2),'b*');
errorbar([3:5],allmeanpix4(3:5),allstdpix4(3:5),'g*');
errorbar([6:8],allmeanpix4(6:8),allstdpix4(6:8),'b*');
xlabel('recording number');
ylabel('mean IBI during recording');
legend('Spontaneous (No Perfusion) (0:00)','Spontaneous (Perfusion (0.5 mL/min)) (0:03)', 'Spontaneous (No Perfusion) (0:13-0:40)', 'Spontaneous (Perfusion (0.5 mL/min,0.75mL/min,1.0mL/min) (0:51-1:08))');

subplot(4,1,3);
recnum=1:8;
Therm=[35 35 36.5 36.5 34.5 34.5 34 33.5];
RealTemp=[36 36 38 38 35.5 35.5 35 34.5];
msize=15;
plot(recnum(1),Therm(1),'g*','MarkerSize',msize); ylim([33 39]);
grid on;
xlim([1 8.5])
hold;
plot(recnum(1),RealTemp(1),'go','MarkerSize',msize);
plot(recnum(2),Therm(2),'b*','MarkerSize',msize);
plot(recnum(2),RealTemp(2),'bo','MarkerSize',msize);
plot(recnum(3:5),Therm(3:5),'g*','MarkerSize',msize);
plot(recnum(3:5),RealTemp(3:5),'go','MarkerSize',msize);
plot(recnum(6:8),Therm(6:8),'b*','MarkerSize',msize);
plot(recnum(6:8),RealTemp(6:8),'bo','MarkerSize',msize);
ylabel('Temperature','FontSize',16);

subplot(4,1,4);
perfrate=[0 1.15 0 0 0 1.15 1.725 2.3];
plot(recnum(1),perfrate(1),'g*');
hold on;
plot(recnum(2),perfrate(2),'b*');
plot(recnum(3:5),perfrate(3:5),'g*');
plot(recnum(6:8),perfrate(6:8),'b*');
grid on;
xlim([1 8.5])
xlabel('Recording Number','FontSize',16); ylabel('Perfusion Rate (mL/min)','FontSize',16);

hold off

figure(5);

recnum=1:8;
Therm=[35 35 36.5 36.5 34.5 34.5 34 33.5];
RealTemp=[36 36 38 38 35.5 35.5 35 34.5];
msize=15;
plot(recnum(1),Therm(1),'g*','MarkerSize',msize); ylim([33 39]);
grid on;
hold;
plot(recnum(1),RealTemp(1),'go','MarkerSize',msize);
plot(recnum(2),Therm(2),'b*','MarkerSize',msize);
plot(recnum(2),RealTemp(2),'bo','MarkerSize',msize);
plot(recnum(3:5),Therm(3:5),'g*','MarkerSize',msize);
plot(recnum(3:5),RealTemp(3:5),'go','MarkerSize',msize);
plot(recnum(6:8),Therm(6:8),'b*','MarkerSize',msize);
plot(recnum(6:8),RealTemp(6:8),'bo','MarkerSize',msize);
legend('Spontaneous (No Perfusion) Thermister Reading (0:00)','Spontaneous (No Perfusion) Real Temperature (0:00)','Spontaneous (Perfusion (0.5 mL/min)) Thermister Reading (0:03)', 'Spontaneous (Perfusion (0.5 mL/min)) Real Temperature (0:03)','Spontaneous (No Perfusion) Thermister Reading (0:13-0:40)', 'Spontaneous (No Perfusion) Real Temperature (0:13-0:40)','Spontaneous (Perfusion (0.5 mL/min,0.75mL/min,1.0mL/min) Thermister Reading (0:51-1:08))','Spontaneous (Perfusion (0.5 mL/min,0.75mL/min,1.0mL/min) Real Temperature (0:51-1:08))');
xlabel('Recording Number'); ylabel('Temperature'); title('May 13 Dish2');

figure(6);
RealTemp=[37 36.7 36.2 36 35 34 33 32 30 29 27];
ThermTemp=[35.5 35.4 35 35 34 33 32 32 30 29 27];
plot(ThermTemp,RealTemp,'g*-');
grid on;
xlabel('Thermister Reading','FontSize',18); ylabel('Calibrated Tempurature','FontSize',18); title('Calibration Curve for Thermister','FontSize',30);