%% Time series for one pixel, simulations

%Assume m holds the signal for one pixel
time=[28.97/66640:28.97/66640:28.97]';
plot(time,squeeze(m(40,40,:)),'LineWidth',2)
xlabel('Time (s.)')
ylabel('Voltage (V)')

%%