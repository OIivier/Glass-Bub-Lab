% Made by MinJu You in 2016

function [SF2] = bandpass(SF1)

% What does it do? 
% 1. Finds the peak frequencies in SF1 through FFT.  
% 2. Filters for a band of peak f through convolution of the data with a filter based on a Hamming Window (fir1).
% Filtered data is stored in SF2. 

% Any time you have an index number instead of frequency, you can convert it to frequency by doing: (Index)*(Fs/nfft2). 

tic

%% Secion 1. Finds the normalized peak frequencies through FFT.  
temp = SF1;

nx = size(temp,1);
ny = size(temp,2);
nframe = size(temp,3);
Fs = 40; % sampling frequency 40Hz 
nfft = nframe; % length of the time domain signal = dt
nfft2 = 2^nextpow2(nfft); 
% nextpow2 = returns the first P such that 2.^P >= abls(N). Ie. if signal length is 2000 = ~2^10.97, then fft duration should be 2048=2^11.
% To get a good resolution in fft, you have to sample at the next power of 2 from the signal length. 

Norm_peakf = zeros(nx,ny);

for ii = 1:nx
    for jj = 1:ny
        x = squeeze(temp(ii,jj,:));
        Y = fft(x, nfft2);
        P0 = Y(1:(nfft2/2));
        [~,Peakf] = max(P0); % Peakf is the index of Y when it is max.
        
        % This section is used to reselect Peakf if they are unreasonable.
        % For example, values that are too small (too fast) or too high
        % (too slow) are removed.
        maxPeakf = 20; % arbitrary value 
        if Peakf >= 1 && Peakf < maxPeakf
           [~, second_max] = max(P0(maxPeakf:nfft2/4)); % nfft2/2 would be normally sufficient, but noisy data can have really high f noise
           Peakf = second_max+(maxPeakf-1); % You have to add (50-1) bc the new index has shifted up by (50-1). 
           
           if Peakf == maxPeakf % This means that most likley there is only noise in this pixel. 
              Peakf = 10000; % 48.8281Hz; arbitrary value given to noisy pixel. 
           end
        end
        Norm_peakf(ii,jj) = (Peakf)/(nfft2/2);
    end
end 
    
%% Section 2: Filters through bandpass and stores filtered data into SF2.  

order = 32; 
SF2_temp = zeros(nx,ny,nframe+order); % ****** What does this last value actually represent? nframe + order! ****** 

for ii = 1:nx
    for jj = 1:ny
        x = squeeze(temp(ii,jj,:));
        
        cutoff = Norm_peakf(ii,jj)*0.05; %0.10
        % Cutoff value has been set arbitrarily. 
        % If you want to compare Cutoff to Hz -> cutoff*(nfft2/2)*(Fs/nfft2)
        empt_pix = 1/(nfft2/2); 
        % Empty pixels that are masked have a Peakf of 1. These pixels are not bandpassed. 

        if cutoff > empt_pix 
            below = cutoff; 
        elseif cutoff < empt_pix 
            below = empt_pix; 
        end

        if Norm_peakf(ii,jj) > below && Norm_peakf(ii,jj) < (1-cutoff) 
            low = Norm_peakf(ii,jj)-cutoff; % cut offs normalized to nyquist frequency so that the value lies between 0 and 1.  
            high = Norm_peakf(ii,jj)+cutoff;
            bandpass_filter = fir1(order,[low high],'bandpass');
            
            % Conv function multiplies the 2 functions to amplify the peaks
            % and dim the troughs and to return the signal into the time 
            % domain.
            con = conv(x,bandpass_filter);             
        else
            % Any frequency that does not lie within the desired range is 
            % not converted back to the frequency domain. 
            con = 0;
        end
        
        SF2_temp(ii,jj,:) = con;
        
    end
end

SF2 = SF2_temp(:,:,1+(order/2):nframe+(order/2)); 
% POTENTIAL FOR ERROR: The final SF2 file that is created has more frames 
% than the original SF1, so this has to be taken into account when making
% the video. To prevent this, we remove the last # frames equal to the 
% value of the order.

toc