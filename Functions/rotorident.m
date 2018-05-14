%S = load('C:\Users\Penisbleu\Documents\MATLAB\Simulationsfolder\Simul2_0-0-15_0-0-0_0-0-4.mat');
%x = S.m(30,30,60000:end);
%plot(squeeze(x));
%% Finding peak frequency at probes, from bandpass MinJu code
function rotorboolean = rotorident(m)
temp = m;

nx = size(temp,1);
ny = size(temp,2);
nframe = round(size(temp,3)*1/10,0);
Fs = 2221; % sampling frequency 40Hz 
nfft = nframe; % length of the time domain signal = dt
nfft2 = 2^nextpow2(nfft); 
% nextpow2 = returns the first P such that 2.^P >= abls(N). Ie. if signal length is 2000 = ~2^10.97, then fft duration should be 2048=2^11.
% To get a good resolution in fft, you have to sample at the next power of 2 from the signal length. 

Norm_peakf = zeros(3,3);

dx = round(nx/11,0);
dy = round(ny/11,0);
startrec = round(size(temp,3)*9/10,0);
for ii = 1:10
    for jj = 1:10
        x = squeeze(temp(dx*ii,dy*jj,startrec:end));
        Y = fft(x, nfft2);
        P0 = Y(1:(nfft2/2));
        [~,Peakf] = max(P0); % Peakf is the index of Y when it is max.
        
        % This section is used to reselect Peakf if they are unreasonable.
        % For example, values that are too small (too fast) or too high
        % (too slow) are removed.
        maxPeakf = 100; % arbitrary value 
        if Peakf >= 1 && Peakf < maxPeakf
           [~, second_max] = max(P0(maxPeakf:nfft2/4)); % nfft2/2 would be normally sufficient, but noisy data can have really high f noise
           Peakf = second_max+(maxPeakf-1); % You have to add (50-1) bc the new index has shifted up by (50-1). 
           
           if Peakf == maxPeakf % This means that most likley there is only noise in this pixel. 
              Peakf = 0; % 48.8281Hz; arbitrary value given to noisy pixel. 
           end
        end
        Norm_peakf(ii,jj) = (Peakf)/(nfft2/2);
    end
end 

%% Deciding if rotor or not

rotorboolean = 1;
error = 0.15;
% % disp(Norm_peakf(4,4));
% plot(squeeze(m(5*dx,4*dy,:)));
% hold on
% % disp(Norm_peakf(6,6));
% plot(squeeze(m(5*dx,3*dy,:)));
%booleantensor = zeros(3,3);
%Norm_peakf(2,2) = 1;
%disp(Norm_peakf);
disp(Norm_peakf);
for j = 1:10
    for k = 1:10
%         if Norm_peakf(j,k)<1/2000
%             Norm_peakf(j,k) = 0;
        if Norm_peakf(j,k) ~= 0
            for v = 1:10
                for w = 1:10
                    if (Norm_peakf(v,w)-Norm_peakf(j,k))/Norm_peakf(j,k) > error && Norm_peakf(v,w) ~= 0
                        rotorboolean = 0;
                        %booleantensor(j,k) = 1;
                    end
                end
            end
        end
    end
end
if isequal(Norm_peakf,zeros(10,10))
    rotorboolean = 0;
end
disp(rotorboolean);
% Heatmap2(Norm_peakf,0.029);
             